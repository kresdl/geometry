package com.kresdl.geometry;

import java.awt.Point;
import java.io.Serializable;

/**
 * 2 element vector
 */
@SuppressWarnings("serial")
public class Vec2 implements Serializable {

	/**
	 * Element
	 */
	public double x, y;

	/**
	 * Constructs a vector with all 0s.
	 */
	public Vec2() {
	}

	/**
	 * Constructs a custom vector.
	 * 
	 * @param x x
	 * @param y y
	 */
	public Vec2(double x, double y) {
		this.x = x;
		this.y = y;
	}

	/**
	 * Constructs a vector equal to the vector passed as argument.
	 * 
	 * @param a vecor
	 */
	public Vec2(Vec2 a) {
		x = a.x;
		y = a.y;
	}

	/**
	 * Constructs a vector equal to the point passed as argument.
	 * 
	 * @param a point
	 */
	public Vec2(Point a) {
		x = a.x;
		y = a.y;
	}

	/**
	 * Returns a float array representation of this vector.
	 *
	 * @return float array of size 2.
	 */
	public float[] toFloats() {
		float[] f = new float[2];
		f[0] = (float) x;
		f[1] = (float) y;
		return f;
	}

	/**
	 * Returns the normalized vector.
	 * 
	 * @param a vector to normalize
	 * @return vector
	 */
	public static Vec2 nrm(Vec2 a) {
		double m = Math.sqrt(a.x * a.x + a.y * a.y);
		return new Vec2(a.x / m, a.y / m);
	}

	/**
	 * Subtracts two vectors and returns the resulting vector.
	 * 
	 * @param a vector a
	 * @param b vector b
	 * @return vector
	 */
	public static Vec2 sub(Vec2 a, Vec2 b) {
		return new Vec2(a.x - b.x, a.y - b.y);
	}

	/**
	 * Adds two vectors and returns the resulting vector.
	 * 
	 * @param a vector a
	 * @param b vector b
	 * @return vector
	 */
	public static Vec2 add(Vec2 a, Vec2 b) {
		return new Vec2(a.x + b.x, a.y + b.y);
	}

	/**
	 * Divides a vector by a scalar and returns the resulting vector.
	 * 
	 * @param a vector
	 * @param d divisor
	 * @return vector
	 */
	public static Vec2 div(Vec2 a, double d) {
		return new Vec2(a.x / d, a.y / d);
	}

	/**
	 * Multiplies a vector by a scalar and returns the resulting vector.
	 * 
	 * @param a vector
	 * @param p factor
	 * @return vector
	 */
	public static Vec2 mul(Vec2 a, double p) {
		return new Vec2(a.x * p, a.y * p);
	}

	/**
	 * Returns the cross product of two vectors.
	 * 
	 * @param a vector a
	 * @param b vector b
	 * @return vector
	 */
	public static double cross(Vec2 a, Vec2 b) {
		return a.x * b.y - a.y * b.x;
	}

	/**
	 * Returns the dot product of two vectors.
	 * 
	 * @param a vector a
	 * @param b vector b
	 * @return vector
	 */
	public static double dot(Vec2 a, Vec2 b) {
		return a.x * b.x + a.y * b.y;
	}

	/**
	 * Sets the vector to the values passed and returns itself.
	 * 
	 * @param x x
	 * @param y y
	 * @return this vector
	 */
	public Vec2 use(double x, double y) {
		this.x = x;
		this.y = y;
		return this;
	}

	/**
	 * Returns the vector's length.
	 * 
	 * @return length
	 */
	public double length() {
		return Math.sqrt(x * x + y * y);
	}
}
