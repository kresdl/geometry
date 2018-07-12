package com.kresdl.geometry;

import java.io.Serializable;

/**
 * A 4x4 matrix
 */
@SuppressWarnings("serial")
public class Matrix implements Serializable {

    /**
     * Element
     */
    public double xx, xy, xz, xw, yx, yy, yz, yw, zx, zy, zz, zw, wx, wy, wz, ww;

    /**
     * Constructs an identity matrix.
     */
    public Matrix() {
        xx = 1;
        yx = 0;
        zx = 0;
        wx = 0;

        xy = 0;
        yy = 1;
        zy = 0;
        wy = 0;

        xz = 0;
        yz = 0;
        zz = 1;
        wz = 0;

        xw = 0;
        yw = 0;
        zw = 0;
        ww = 1;
    }

    /**
     * Constructs a custom matrix.
     *
     * @param x xx, yx, zx
     * @param y xy, yy, zy
     * @param z xz, yz, zz
     * @param tx wx
     * @param ty wy
     * @param tz wz
     */
    public Matrix(Vec x, Vec y, Vec z, double tx, double ty, double tz) {
        xx = x.x;
        yx = x.y;
        zx = x.z;
        wx = tx;

        xy = y.x;
        yy = y.y;
        zy = y.z;
        wy = ty;

        xz = z.x;
        yz = z.y;
        zz = z.z;
        wz = tz;

        xw = 0;
        yw = 0;
        zw = 0;
        ww = 1;
    }

    /**
     * Constructs a custom matrix.
     *
     * @param x xx, yx, zx
     * @param y xy, yy, zy
     * @param z xz, yz, zz
     */
    public Matrix(Vec x, Vec y, Vec z) {
        xx = x.x;
        yx = x.y;
        zx = x.z;
        wx = 0;

        xy = y.x;
        yy = y.y;
        zy = y.z;
        wy = 0;

        xz = z.x;
        yz = z.y;
        zz = z.z;
        wz = 0;

        xw = 0;
        yw = 0;
        zw = 0;
        ww = 1;
    }

    /**
     * Constructs a custom matrix.
     *
     * @param x xx, yx, zx
     * @param y xy, yy, zy
     * @param z xz, yz, zz
     * @param p wx, wy, wz
     */
    public Matrix(Vec x, Vec y, Vec z, Vec p) {
        xx = x.x;
        yx = x.y;
        zx = x.z;
        wx = p.x;

        xy = y.x;
        yy = y.y;
        zy = y.z;
        wy = p.y;

        xz = z.x;
        yz = z.y;
        zz = z.z;
        wz = p.z;

        xw = 0;
        yw = 0;
        zw = 0;
        ww = 1;
    }

    /**
     * Constructs a matrix identical to the matrix passed as argument.
     *
     * @param a matrix
     */
    public Matrix(Matrix a) {
        xx = a.xx;
        yx = a.yx;
        zx = a.zx;
        wx = a.wx;

        xy = a.xy;
        yy = a.yy;
        zy = a.zy;
        wy = a.wy;

        xz = a.xz;
        yz = a.yz;
        zz = a.zz;
        wz = a.wz;

        xw = a.xw;
        yw = a.yw;
        zw = a.zw;
        ww = a.ww;
    }

    /**
     * Returns a representation of the matrix as an array of floats.
     *
     * @return a float array with size 16
     */
    public float[] toFloats() {
        float[] f = new float[4 * 4];
        f[0] = (float) xx;
        f[1] = (float) yx;
        f[2] = (float) zx;
        f[3] = (float) wx;

        f[4] = (float) xy;
        f[5] = (float) yy;
        f[6] = (float) zy;
        f[7] = (float) wy;

        f[8] = (float) xz;
        f[9] = (float) yz;
        f[10] = (float) zz;
        f[11] = (float) wz;

        f[12] = (float) xw;
        f[13] = (float) yw;
        f[14] = (float) zw;
        f[15] = (float) ww;
        return f;
    }

    /**
     * Returns a matrix based on rotations around X- Y- and Z axis.
     *
     * @param x rotation angle around X axis, in radians.
     * @param y rotation angle around Y axis, in radians.
     * @param z rotation angle around Z axis, in radians.
     * @return matrix
     *
     */
    public static Matrix r(double x, double y, double z) {
        Matrix m1 = rx(x);
        Matrix m2 = ry(y);
        Matrix m3 = rz(z);
        return mul(mul(m1, m2), m3);
    }

    /**
     * Returns a matrix based on rotation around an arbitrary axis by a specific
     * angle.
     *
     * @param quat a vector consisting of axis.z, axis.y, axiz.w and rotation
     * angle in radians
     * @return matrix
     */
    public static Matrix r(Vec quat) {
        Vec x = Vec.nrm(Vec.cross(quat, new Vec(0, 1, 0)));
        Vec y = Vec.cross(x, quat);
        Matrix v = new Matrix(x, y, quat);
        return mul(mul(v, rz(quat.w)), inv(v));
    }

    /**
     * Returns a matrix based on rotation around x-axis by an angle of r
     * radians.
     *
     * @param r angle in radians
     * @return matrix
     */
    public static Matrix rx(double r) {
        Matrix m = new Matrix();
        double sin = Math.sin(r);
        double cos = Math.cos(r);
        m.yy = cos;
        m.zy = -sin;

        m.yz = sin;
        m.zz = cos;
        return m;
    }

    /**
     * Returns a matrix based on rotation around y-axis by an angle of r
     * radians.
     *
     * @param r angle in radians
     * @return matrix
     */
    public static Matrix ry(double r) {
        Matrix m = new Matrix();
        double sin = Math.sin(r);
        double cos = Math.cos(r);
        m.xx = cos;
        m.zx = -sin;

        m.xz = sin;
        m.zz = cos;
        return m;
    }

    /**
     * Returns a matrix based on rotation around z-axis by an angle of r
     * radians.
     *
     * @param r angle in radians
     * @return matrix
     */
    public static Matrix rz(double r) {
        Matrix m = new Matrix();
        double sin = Math.sin(r);
        double cos = Math.cos(r);
        m.xx = cos;
        m.yx = -sin;

        m.xy = sin;
        m.yy = cos;
        return m;
    }

    /**
     * Returns a scaling matrix.
     *
     * @param s scale
     * @return matrix
     */
    public static Matrix s(double s) {
        Matrix m = new Matrix();
        m.xx = s;
        m.yy = s;
        m.zz = s;
        return m;
    }

    /**
     * Returns a scaling matrix.
     *
     * @param x scale in x-axis
     * @param y scale in y-axis
     * @param z scale in z-axis
     * @return matrix
     */
    public static Matrix s(double x, double y, double z) {
        Matrix m = new Matrix();
        m.xx = x;
        m.yy = y;
        m.zz = z;
        return m;
    }

    /**
     * Returns a translation matrix.
     *
     * @param x translation in x-axis
     * @param y translation in y-axis
     * @param z translation in z-axis
     * @return matrix
     */
    public static Matrix t(double x, double y, double z) {
        Matrix m = new Matrix();
        m.wx = x;
        m.wy = y;
        m.wz = z;
        return m;
    }

    /**
     * Returns a translation matrix.
     *
     * @param p translation vector
     * @return matrix
     */
    public static Matrix t(Vec p) {
        Matrix m = new Matrix();
        m.wx = p.x;
        m.wy = p.y;
        m.wz = p.z;
        return m;
    }

    /**
     * Returns a matrix based on rotation around an arbitrary axis by a specific
     * angle.
     *
     * @param a axis point A
     * @param b axis point B
     * @param angle rotation angle in radians.
     * @return matrix
     */
    public static Matrix arb(Vec a, Vec b, double angle) {
        Matrix v = view(a, b);
        return mul(mul(v, rz(angle)), inv(v));
    }

    /**
     * Returns a view matrix.
     *
     * @param pos view position
     * @param target view target
     * @return matrix
     */
    public static Matrix view(Vec pos, Vec target) {
        Vec z = Vec.nrm(Vec.sub(target, pos));
        Vec x = Vec.nrm(Vec.cross(z, new Vec(0, 1, 0)));
        Vec y = Vec.cross(x, z);
        Matrix m = new Matrix();
        m.xx = x.x;
        m.yx = x.y;
        m.zx = x.z;
        m.wx = -Vec.dot(pos, x);

        m.xy = y.x;
        m.yy = y.y;
        m.zy = y.z;
        m.wy = -Vec.dot(pos, y);

        m.xz = z.x;
        m.yz = z.y;
        m.zz = z.z;
        m.wz = -Vec.dot(pos, z);
        return m;
    }

    /**
     * Returns a view matrix adapted to OpenGL space.
     *
     * @param pos view position
     * @param target view target
     * @return matrix
     */
    public static Matrix glView(Vec pos, Vec target) {
        Vec z = Vec.nrm(Vec.sub(pos, target));
        Vec x = Vec.nrm(Vec.cross(new Vec(0, 1, 0), z));
        Vec y = Vec.cross(z, x);
        Matrix m = new Matrix();
        m.xx = x.x;
        m.yx = x.y;
        m.zx = x.z;
        m.wx = -Vec.dot(pos, x);

        m.xy = y.x;
        m.yy = y.y;
        m.zy = y.z;
        m.wy = -Vec.dot(pos, y);

        m.xz = z.x;
        m.yz = z.y;
        m.zz = z.z;
        m.wz = -Vec.dot(pos, z);
        return m;
    }

    /**
     * Returns a projection matrix.
     *
     * @param zNear z near
     * @param zFar z far
     * @param fov angle of field-of-view in radians.
     * @param aspectRatio aspect ratio (viewport width / wiewport height)
     * @return matrix
     */
    public static Matrix proj(double zNear, double zFar, double fov, double aspectRatio) {
        Matrix m = new Matrix();
        double py = 1.0d / Math.tan(fov / 2);
        double px = py / aspectRatio;
        double pz = zFar / (zFar - zNear);
        m.xx = px;

        m.yy = -py;

        m.zz = pz;
        m.wz = -pz * zNear;

        m.zw = 1;
        m.ww = 0;
        return m;
    }

    /**
     * Returns a projection matrix adapted to OpenGL space.
     *
     * @param zNear z near
     * @param zFar z far
     * @param fov angle of field-of-view in radians.
     * @param aspectRatio aspect ratio (viewport width / wiewport height)
     * @return matrix
     */
    public static Matrix glProj(double zNear, double zFar, double fov, double aspectRatio) {
        Matrix m = new Matrix();
        double py = 1.0d / Math.tan(fov / 2);
        double px = py / aspectRatio;
        double pz = (zFar + zNear) / (zNear - zFar);
        m.xx = px;

        m.yy = py;

        m.zz = pz;
        m.wz = 2 * zNear * zFar / (zNear - zFar);

        m.zw = -1;
        m.ww = 0;
        return m;
    }

    /**
     * Returns the concatenated matrix from two matrices.
     *
     * @param a matrix a
     * @param b matrix b
     * @return matrix
     */
    public static Matrix mul(Matrix a, Matrix b) {
        Matrix m = new Matrix();
        m.xx = (b.xx * a.xx) + (b.yx * a.xy) + (b.zx * a.xz) + (b.wx * a.xw);
        m.yx = (b.xx * a.yx) + (b.yx * a.yy) + (b.zx * a.yz) + (b.wx * a.yw);
        m.zx = (b.xx * a.zx) + (b.yx * a.zy) + (b.zx * a.zz) + (b.wx * a.zw);
        m.wx = (b.xx * a.wx) + (b.yx * a.wy) + (b.zx * a.wz) + (b.wx * a.ww);

        m.xy = (b.xy * a.xx) + (b.yy * a.xy) + (b.zy * a.xz) + (b.wy * a.xw);
        m.yy = (b.xy * a.yx) + (b.yy * a.yy) + (b.zy * a.yz) + (b.wy * a.yw);
        m.zy = (b.xy * a.zx) + (b.yy * a.zy) + (b.zy * a.zz) + (b.wy * a.zw);
        m.wy = (b.xy * a.wx) + (b.yy * a.wy) + (b.zy * a.wz) + (b.wy * a.ww);

        m.xz = (b.xz * a.xx) + (b.yz * a.xy) + (b.zz * a.xz) + (b.wz * a.xw);
        m.yz = (b.xz * a.yx) + (b.yz * a.yy) + (b.zz * a.yz) + (b.wz * a.yw);
        m.zz = (b.xz * a.zx) + (b.yz * a.zy) + (b.zz * a.zz) + (b.wz * a.zw);
        m.wz = (b.xz * a.wx) + (b.yz * a.wy) + (b.zz * a.wz) + (b.wz * a.ww);

        m.xw = (b.xw * a.xx) + (b.yw * a.xy) + (b.zw * a.xz) + (b.ww * a.xw);
        m.yw = (b.xw * a.yx) + (b.yw * a.yy) + (b.zw * a.yz) + (b.ww * a.yw);
        m.zw = (b.xw * a.zx) + (b.yw * a.zy) + (b.zw * a.zz) + (b.ww * a.zw);
        m.ww = (b.xw * a.wx) + (b.yw * a.wy) + (b.zw * a.wz) + (b.ww * a.ww);
        return m;
    }

    /**
     * Returns the matrix inverse.
     *
     * @param m matrix to invert
     * @return matrix
     */
    public static Matrix inv(Matrix m) {
        Matrix a = new Matrix();
        Vec p = new Vec(m.wx, m.wy, m.wz);
        a.xx = m.xx;
        a.yx = m.xy;
        a.zx = m.xz;
        a.wx = -Vec.dot(p, new Vec(m.xx, m.xy, m.xz));

        a.xy = m.yx;
        a.yy = m.yy;
        a.zy = m.yz;
        a.wy = -Vec.dot(p, new Vec(m.yx, m.yy, m.yz));

        a.xz = m.zx;
        a.yz = m.zy;
        a.zz = m.zz;
        a.wz = -Vec.dot(p, new Vec(m.zx, m.zy, m.zz));
        return a;
    }

    /**
     * Returns the matrix transpose.
     *
     * @param m matrix to transpose
     * @return matrix
     */
    public static Matrix tp(Matrix m) {
        Matrix t = new Matrix();
        t.xx = m.xx;
        t.yx = m.xy;
        t.zx = m.xz;
        t.wx = m.xw;

        t.xy = m.yx;
        t.yy = m.yy;
        t.zy = m.yz;
        t.wy = m.yw;

        t.xz = m.zx;
        t.yz = m.zy;
        t.zz = m.zz;
        t.wz = m.zw;

        t.xw = m.wx;
        t.yw = m.wy;
        t.zw = m.wz;
        t.ww = m.ww;
        return t;
    }
}
