package ingo;

import robocode.*;
import java.util.Hashtable;
import java.io.*;
import java.awt.Color;

public class Ingoo extends AdvancedRobot {

	public static final boolean rnew = false;
	static int[] deadbots;
	Surf move;
	Radar radar;
	TankTube gun;

	public void run() {
		init();
		int i = 0;
		while (true) {
			if (!rnew)
				try {
					move.onTick();
				} catch (Exception ex) {
					ex.printStackTrace();
				}
			try {
				radar.onTick();
			} catch (Exception ex) {
				ex.printStackTrace();
			}
			try {
				gun.onTick();
			} catch (Exception ex) {
				ex.printStackTrace();
			}
			if (i++ >= 10)
				try {
					i = 0;
					Thread.sleep(0, 1);
				} catch (Exception e) {
				}
			execute();
		}
	}

	public void init() {
		getDataDirectory();
		if (deadbots == null)
			deadbots = new int[getOthers() + 1];
		setAdjustRadarForGunTurn(true);
		setAdjustRadarForRobotTurn(true);
		setAdjustGunForRobotTurn(true);
		setColors(Color.BLACK, Color.PINK, Color.BLACK);

		if (!rnew) {
			try {
				move = Surf.getInstance(this);
			} catch (Exception ex) {
				ex.printStackTrace();
			}
		}
		try {
			radar = new Radar(this);
		} catch (Exception ex) {
			ex.printStackTrace();
		}
		try {
			gun = new TankTube(this);
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}

	public void onSkippedTurn(SkippedTurnEvent e) {

	}

	public void onScannedRobot(ScannedRobotEvent e) {
		if (!rnew)
			try {
				move.onScannedRobot(e);
			} catch (Exception ex) {
				ex.printStackTrace();
			}
		try {
			radar.onScannedRobot(e);
		} catch (Exception ex) {
			ex.printStackTrace();
		}
		try {
			gun.onScannedRobot(e);
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}

	public void onRobotDeath(RobotDeathEvent e) {
		if (!rnew)
			try {
				move.onRobotDeath(e);
			} catch (Exception ex) {
				ex.printStackTrace();
			}
		try {
			radar.onRobotDeath(e);
		} catch (Exception ex) {
			ex.printStackTrace();
		}
		try {
			gun.onRobotDeath(e);
		} catch (Exception ex) {
			ex.printStackTrace();
		}

	}

	public void onHitByBullet(HitByBulletEvent e) {
		if (!rnew)
			try {
				move.onHitByBullet(e);
			} catch (Exception ex) {
				ex.printStackTrace();
			}
	}

	public void onBulletHitBullet(BulletHitBulletEvent e) {
		if (!rnew)
			try {
				move.onBulletHitBullet(e);
			} catch (Exception ex) {
				ex.printStackTrace();
			}
	}

	public void onBulletHit(BulletHitEvent e) {
		if (!rnew)
			try {
				move.onBulletHit(e);
			} catch (Exception ex) {
				ex.printStackTrace();
			}
	}

	public void onPaint(java.awt.Graphics2D g) {
		if (!rnew)
			move.onPaint(g);
		gun.onPaint(g);
	}

	public void onDeath(DeathEvent e) {
		try {
			gun.onDeath(e);
			deadbots[getOthers()]++;
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}

	public void onWin(WinEvent e) {
		try {
			gun.onWin(e);
			deadbots[0]++;
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}

	public Hashtable<String, Double> predictAllBulletPowers() {
		return move.predictAllBulletPowers();
	}

	public void bulletFired(Bullet b) {
		if (!rnew)
			try {
				move.bulletFired(b);
			} catch (Exception ex) {
				ex.printStackTrace();
			}
	}

}