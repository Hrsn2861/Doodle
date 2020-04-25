package cn.hzw.doodle.STFT;

import android.util.Log;

public class HammingWindow {
    static double coef = 0.46;
    public HammingWindow(double coef) { this.coef = coef; }
    public HammingWindow() { }

    public static float[] calHamming(float[] x) {
        String s = "";
        for(int i=0;i<x.length;i++) {
            s += i+": "+((1-coef) - coef * Math.cos(2 * Math.PI * i/(x.length-1)))+", ";
            x[i] *= ((1-coef) - coef * Math.cos(2 * Math.PI * i/(x.length-1)));
        }
        return x;
    }
}
