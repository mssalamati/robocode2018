package ingo;
import java.awt.geom.*;
import java.util.*;

import robocode.util.*;

public class Wave implements HistoryLog.LogListener{

    double waveWeight;

    String firedBy;
    DangerBot firer;

    Point2D.Double fireLocation;

    Hashtable<String,DangerBot> snapshot = new Hashtable<String,DangerBot>();

    long fireTime;

    double bulletPower, bulletVelocity, bulletDamage;
    double firerEnergy;

    double[] bins;
    double[] botShadowBins;

    boolean gunHeatWave = false;
    boolean surfable = true;

    boolean needsSnapshotRebuild = true;
    boolean needsDangerRecalc = true;

    void updateIfNecessary(HistoryLog log, Point2D.Double myLocation){
        if(needsSnapshotRebuild)
            buildSnapshot(log);
        if(needsDangerRecalc)
            calcDangers(myLocation);
    }
    public void newData(String botName){
        needsSnapshotRebuild = true;
    }

    void buildSnapshot(HistoryLog log){
        snapshot.clear();
        Point2D.Double[] corners = new Point2D.Double[]{new Point2D.Double(0,0),
                                                        new Point2D.Double(0,Surf.MAX_Y),
                                                        new Point2D.Double(Surf.MAX_X,0),
                                                        new Point2D.Double(Surf.MAX_X,Surf.MAX_Y)};
 
        long dataTime = fireTime;
        HashMap<String, HistoryLog.InterpolatedLogEntry> botData = log.getAllNearest(dataTime);

        HistoryLog.InterpolatedLogEntry firerAtFire = botData.get(firedBy);
        if (firerAtFire == null)
            return;

        fireLocation = firerAtFire.location;
        firerEnergy = firerAtFire.energy;

        Set<String> botNames = botData.keySet();

        for(String botName : botNames){
            HistoryLog.InterpolatedLogEntry ei = botData.get(botName);
            if(botName.equals(firedBy) || ei == null)
                continue;

            if(ei.waitingOnData)
                log.updateOnNewData(botName, this);

            DangerBot cp = new DangerBot();
            cp.location = ei.location;
            cp.energy = ei.energy;
            cp.name = botName;
            cp.heading = ei.heading;
            cp.velocity = ei.velocity;

            double targetBearing = Surf.absoluteBearing(fireLocation,cp.location);
            cp.latVel = cp.velocity*HandBrake.sin(cp.heading - targetBearing);
            cp.advVel = cp.velocity*HandBrake.cos(cp.heading - targetBearing);
            cp.distToE = cp.location.distance(fireLocation);
            double distToNearestSq = cp.distToE * cp.distToE;

            for(String botNameNearest : botNames){
                HistoryLog.InterpolatedLogEntry ein = botData.get(botNameNearest);
                if(botNameNearest.equals(cp.name))
                    continue;
                double distSq = ein.location.distanceSq(cp.location);
                if(distSq < distToNearestSq){
                    distToNearestSq = distSq;
                }
            }
            cp.distToNearest = Math.sqrt(distToNearestSq);

            cp.distToWall = Math.min(
                            Math.min(cp.location.x - 18,cp.location.y - 18),
                            Math.min(Surf.MAX_X - 18 - cp.location.x, Surf.MAX_Y - 18 - cp.location.y));

            cp.enemiesAlive = botNames.size();

            double distToCornerSq = Double.POSITIVE_INFINITY;
            for(int i = 0; i < 4; i++)
                distToCornerSq = Math.min(cp.location.distanceSq(corners[i]),distToCornerSq);
            cp.distToCorner = Math.sqrt(distToCornerSq);

            //acceleration
            HistoryLog.InterpolatedLogEntry back1 = log.getInterpolatedNearest(botName, dataTime-1);
            cp.accel = Math.abs(ei.velocity) - Math.abs(back1.velocity);

            //dist-last-10
            HistoryLog.InterpolatedLogEntry back10 = log.getInterpolatedNearest(botName, dataTime-10);
            cp.distLast10 = back10.location.distance(ei.location);
            double dir = Math.signum(cp.velocity);
            for(int i = 1; i < 30; i++){
                HistoryLog.InterpolatedLogEntry backi = log.getInterpolatedNearest(botName, dataTime-i);
                if (Math.signum(backi.velocity) != dir)
                    break;
                cp.timeSinceReverse++;
            }

            snapshot.put(cp.name,cp);
        }
        needsDangerRecalc = true;
        needsSnapshotRebuild = false;
    }
    double sqr(double s){
        return s*s;
    }

    public void calcDangers(Point2D.Double myLocation){
        surfable = false;
        Enumeration<DangerBot> en = snapshot.elements();
        while(en.hasMoreElements()){
            DangerBot target = en.nextElement();

            if(target.name.equals(firedBy))
                continue;

            double targetBearing = Surf.absoluteBearing(fireLocation,target.location);
            double latVel = target.velocity*HandBrake.sin(target.heading - targetBearing);
            double distToE = target.location.distance(fireLocation);
            double meBearing = Surf.absoluteBearing(fireLocation,myLocation);

            KDTree.WeightedManhattan<Scanner> tree = firer.targets.get(target.name);
            int k = 0;
            if(tree == null)
                if(firer.defaultAim != null){
                    tree = firer.defaultAim;
                    k = 1;
                }
                else{
                    tree = Surf.GF_0_tree;
                    k = 1;
                }

            double MEA = Surf.maxEscapeAngle(bulletVelocity);

            if(Math.abs(Utils.normalRelativeAngle(targetBearing - meBearing)) > 3*MEA)
                continue;

            double GFcorrection = Math.signum(latVel)*MEA;

            double botWeight = 1/(distToE*distToE);
            double[] bins = new double[360];

            do{
                List<KDTree.SearchResult<Scanner>> cluster = tree.nearestNeighbours(
                            target.targetDescriptor(),
                            Math.min(10,tree.size())
                            );

                Iterator<KDTree.SearchResult<Scanner>> it = cluster.iterator();
                double weight = botWeight*Math.exp(-k);
                double clusterDistance = 0;
                while(it.hasNext()){
                    KDTree.SearchResult<Scanner> v = it.next();
                    clusterDistance += v.distance;
                }
                clusterDistance /= cluster.size();
                it = cluster.iterator();
                while(it.hasNext()){
                    KDTree.SearchResult<Scanner> v = it.next();
                    double fireAngle = Utils.normalAbsoluteAngle
                                    (v.payload.GF*GFcorrection + targetBearing);

                    if(Math.abs(Utils.normalRelativeAngle(fireAngle - meBearing)) < 2*MEA){
                        Surf.smoothAround(bins,((int)(fireAngle*(180/Math.PI)))%360,
                                    18,v.payload.weight*weight*HandBrake.exp(-0.5*sqr(v.distance/clusterDistance)));
                        surfable = true;
                    }
                }
                k++;

                tree = firer.defaultAim;

            }while(k == 1);


            Surf.areaNormalize(bins);
            for(int i = 0; i < bins.length; i++)
                this.bins[i] += botWeight*bins[i];
        }

        Surf.areaNormalize(this.bins);
        needsDangerRecalc = false;
    }

    void logBulletForShadows(Point2D.Double currentbP, double heading, double velocity, long time)
    {
        if (!surfable)
            return;

        if (sqr(bulletVelocity*(time - fireTime)) > currentbP.distanceSq(fireLocation))
            return;

        double t = 0;
        while(t < 91){
            t++;
            Point2D.Double bP = Surf.project(currentbP,heading,velocity*(t+1));
            Point2D.Double lastbP = Surf.project(currentbP,heading,velocity*t);
            double waveRadius = bulletVelocity * (time - fireTime + t);
            double lastWaveRadius = waveRadius - bulletVelocity;

            double en = bP.distance(fireLocation);
            double len = lastbP.distance(fireLocation);

            if (en < lastWaveRadius || len < en)
                break;

            if (en < waveRadius && len > lastWaveRadius && en < len){

                Point2D.Double p1, p2;

                if(en >= lastWaveRadius)// 1 & 4
                    p1 = bP;
                else{// 2 & 3
                    PreciseWave wv = new PreciseWave();
                    wv.fireLocation = fireLocation;
                    wv.distanceTraveled = lastWaveRadius;
                    p1 = PreciseUtils.intersection(lastbP,bP,wv);
                }
                if(len > waveRadius){//1 & 3
                    PreciseWave wv = new PreciseWave();
                    wv.fireLocation = fireLocation;
                    wv.distanceTraveled = waveRadius;
                    p2 = PreciseUtils.intersection(lastbP,bP,wv);
                }
                else // 2 & 4
                    p2 = lastbP;

                double a1 = Surf.absoluteBearing(fireLocation, p1);
                double a2 = Surf.absoluteBearing(fireLocation, p2);
                double aDiff = Utils.normalRelativeAngle(a2 - a1);
                double angle = Utils.normalAbsoluteAngle(a1 + aDiff/2);
                double width = Math.abs(aDiff);
                logShadow(angle, width);
            }
        }
    }
    void logShadow(double angle, double width){

        double lowIndexD = bins.length*Utils.normalAbsoluteAngle(angle - 0.5*width)*(0.5/Math.PI);
        double highIndexD = bins.length*Utils.normalAbsoluteAngle(angle + 0.5*width)*(0.5/Math.PI);
        int lowIndex = (int)Math.ceil(lowIndexD)%bins.length;
        int highIndex = (int)Math.floor(highIndexD)%bins.length;
        if(lowIndexD <= highIndexD)
            for(int i = lowIndex; i <= highIndex; i++)
                botShadowBins[i] = 0;
        else {
            for(int i = lowIndex; i < botShadowBins.length; i++)
                botShadowBins[i] = 0;
            for(int i = 0; i <= highIndex; i++)
                botShadowBins[i] = 0;
        }
    }
    void checkShadows(Point2D.Double bulletLocation){
        double bulletAngle = Utils.normalAbsoluteAngle(Surf.absoluteBearing(fireLocation, bulletLocation));
        int index = (int)(botShadowBins.length*bulletAngle*0.5/Math.PI)%botShadowBins.length;
        if (botShadowBins[index] == 0.0)
            System.out.println("Hit in bullet shadow!");
    
    }

}