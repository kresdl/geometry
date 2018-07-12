package com.kresdl.geometry;

import java.io.Serializable;

/**
 * 4 element vector
 */
@SuppressWarnings("serial")
public class Vec implements Serializable {

    /**
     * Element
     */
    public double x, y, z, w = 0;

    /**
     * Adds two vectors and returns the resulting vector.
     *
     * @param a vector a
     * @param b vector b
     * @return vector
     */
    public static Vec add(Vec a, Vec b) {
        return new Vec(a.x + b.x, a.y + b.y, a.z + b.z);
    }

    /**
     * Returns the cross product of two vectors.
     *
     * @param a vector a
     * @param b vector b
     * @return vector
     */
    public static Vec cross(Vec a, Vec b) {
        double x = a.y * b.z - a.z * b.y;
        double y = a.z * b.x - a.x * b.z;
        double z = a.x * b.y - a.y * b.x;
        return new Vec(x, y, z);
    }

    /**
     * Divides a vector by a scalar and returns the resulting vector.
     *
     * @param a vector
     * @param b divisor
     * @return vector
     */
    public static Vec div(Vec a, double b) {
        return new Vec(a.x / b, a.y / b, a.z / b);
    }

    /**
     * Returns the dot product of two vectors.
     *
     * @param a vector a
     * @param b vector b
     * @return vector
     */
    public static double dot(Vec a, Vec b) {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }

    /**
     * Interpolates two vectors and returns the resulting vector.
     *
     * @param a vector a
     * @param b vector b
     * @param c parameter between 0.0 and 1.0
     * @return vector
     */
    public static Vec lerp(Vec a, Vec b, double c) {
        return new Vec(a.x + c * (b.x - a.x), a.y + c * (b.y - a.y), a.z + c * (b.z - a.z));
    }

    /**
     * Multiplies a vector by a scalar and returns the resulting vector.
     *
     * @param a vector
     * @param b factor
     * @return vector
     */
    public static Vec mul(Vec a, double b) {
        return new Vec(a.x * b, a.y * b, a.z * b);
    }

    /**
     * Returns the normalized vector.
     *
     * @param a vector to normalize
     * @return vector
     */
    public static Vec nrm(Vec a) {
        return div(a, a.length());
    }

    /**
     * Returns the negated vector.
     *
     * @param a vector to negate
     * @return vector
     */
    public static Vec neg(Vec a) {
        return new Vec(-a.x, -a.y, -a.z);
    }

    /**
     * Subtracts two vectors and returns the resulting vector.
     *
     * @param a vector a
     * @param b vector b
     * @return vector
     */
    public static Vec sub(Vec a, Vec b) {
        return new Vec(a.x - b.x, a.y - b.y, a.z - b.z);
    }

    /**
     * Transforms a vector, treating it like a coordinate and returns the
     * result.
     *
     * @param coord coordinate
     * @param m matrix
     * @return vector
     */
    public static Vec tc(Vec coord, Matrix m) {
        Vec a = new Vec();
        a.x = coord.x * m.xx + coord.y * m.yx + coord.z * m.zx + m.wx;
        a.y = coord.x * m.xy + coord.y * m.yy + coord.z * m.zy + m.wy;
        a.z = coord.x * m.xz + coord.y * m.yz + coord.z * m.zz + m.wz;
        a.w = coord.x * m.xw + coord.y * m.yw + coord.z * m.zw + m.ww;
        return a;
    }

    /**
     * Transforms a vector, treating it like a normal and returns the result.
     *
     * @param nrm normal
     * @param m matrix
     * @return vector
     */
    public static Vec tn(Vec nrm, Matrix m) {
        Vec a = new Vec();
        a.x = nrm.x * m.xx + nrm.y * m.yx + nrm.z * m.zx;
        a.y = nrm.x * m.xy + nrm.y * m.yy + nrm.z * m.zy;
        a.z = nrm.x * m.xz + nrm.y * m.yz + nrm.z * m.zz;
        a.w = 0;
        return a;
    }

    /**
     * Constructs a vector with all 0s.
     */
    public Vec() {
    }

    /**
     * Constructs a custom vector.
     *
     * @param x x
     * @param y y
     * @param z z
     */
    public Vec(double x, double y, double z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }

    /**
     * Constructs a custom vector.
     *
     * @param x x
     * @param y y
     * @param z z
     * @param w w
     */
    public Vec(double x, double y, double z, double w) {
        this.x = x;
        this.y = y;
        this.z = z;
        this.w = w;
    }

    /**
     * Constructs a vector equal to the vector passed as argument.
     *
     * @param a vecor
     */
    public Vec(Vec a) {
        x = a.x;
        y = a.y;
        z = a.z;
        w = a.w;
    }

    /**
     * Returns the vector's length.
     *
     * @return length
     */
    public double length() {
        return Math.sqrt(x * x + y * y + z * z);
    }

    /**
     * Returns a float array representation of the vector.
     *
     * @return float array of size 3.
     */
    public float[] toFloats() {
        float[] f = new float[3];
        f[0] = (float) x;
        f[1] = (float) y;
        f[2] = (float) z;
        return f;
    }

    /**
     * Sets the vector to the values passed and returns itself.
     *
     * @param x x
     * @param y y
     * @param z z
     * @return this vector
     */
    public Vec use(double x, double y, double z) {
        this.x = x;
        this.y = y;
        this.z = z;
        return this;
    }
}
