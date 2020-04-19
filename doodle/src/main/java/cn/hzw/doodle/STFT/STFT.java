package cn.hzw.doodle.STFT;

import android.util.Log;

import java.io.IOException;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;

/**
 *
 * @author Paco
 *
 */
public class STFT {

    float framedSignal[][];         // 切分之后的framedSignal
    ArrayList<ArrayList<Complex>> magnitude;            // rfft 之后得到的幅度值
    double[][] xdb;                   // 计算 amplitude_to_db 的结果

    void addWinfun() {
        for(int i=0;i<framedSignal.length;i++) {
            framedSignal[i] = HammingWindow.calHamming(framedSignal[i]);
        }
    }

    public double[][] performStft(float[] signal, int sampleRate, double frameSize, double frameStride, int NFFT, boolean normalized) {
        int frameLength = (int)Math.round(frameSize *  sampleRate);
        int frameStep = (int)Math.round(frameStride * sampleRate);
        int signalLength = signal.length;
        int numFrames = 1 + (int)(Math.ceil((float)(Math.abs(signalLength - frameLength)) / frameStep));

        int padLength = (int)(numFrames * frameStep + frameLength);
        float[] padSignal = new float[padLength];
        for(int i=0;i<signalLength;i++) { padSignal[i] = signal[i]; }           // pad zeros to the signal

        framedSignal = new float[numFrames][frameLength];
        magnitude = new ArrayList<ArrayList<Complex>>();

        for(int i=0;i<numFrames;i++) {
            for(int j=0;j<frameLength;j++) {
                framedSignal[i][j] = padSignal[frameStep * numFrames + j];
            }
        }

        addWinfun();        // ad hamming window to the signal

        for(int i=0;i<numFrames;i++) {
            magnitude.add(calculateFFT(framedSignal[i], frameLength, NFFT));
        }

        assert magnitude.size() == numFrames;
        assert magnitude.get(0).size() == (NFFT / 2 - 1);

        return calculateFeature3(normalized);
    }

    double AverageXdb() {
        double ret = 0;
        for(int i=0;i<xdb.length;i++) {
            for(int j=0;j<xdb[i].length;j++) {
                ret += xdb[i][j];
            }
        }
        return ret;
    }

    double StdVariantXdb() {
        int len = xdb.length * xdb[0].length;
        double avg = AverageXdb();
        double dVar = 0;
        for(int i=0;i<xdb.length;i++) {
            for(int j=0;j<xdb[i].length;j++) {
                double xi = xdb[i][j];
                dVar += (xi - avg) * (xi - avg);
            }
        }
        return Math.sqrt(dVar/len);
    }

    public double[][] calculateFeature3(boolean normalized) {
        amplitude_to_db();      // Matrix
        double avg = AverageXdb();
        double std = StdVariantXdb();
        if(normalized) {
            for(int i=0;i<xdb.length;i++) {
                for (int j=0;j<xdb[i].length;j++) {
                    double xi = (xdb[i][j] - avg) / std;
                    xdb[i][j] = xi;
                }
            }
        }
        double[][] ret = new double[xdb[0].length][xdb.length];
        // transpose
        for(int i=0;i<xdb[0].length;i++) {
            for(int j=0;j<xdb.length;j++) {
                ret[i][j] = xdb[j][i];
            }
        }

        return ret;
    }

    void amplitude_to_db() {
        double amin = 1E-10;
        double topDb = 80.0;

        xdb = new double[magnitude.size()][magnitude.get(0).size()];

        double[][] squareMagnitude = new double[magnitude.size()][magnitude.get(0).size()];

        for(int i=0;i<magnitude.size();i++) {
            for(int j=0;j<magnitude.get(0).size();j++) {
                squareMagnitude[i][j] = magnitude.get(i).get(j).square();
            }
        }

        double[][] logSpec = new double[magnitude.size()][magnitude.get(0).size()];
        double maxSpec = 0;
        for(int i=0;i<magnitude.size();i++) {
            for(int j=0;j<magnitude.get(0).size();j++) {
                double d = 10 * Math.log10(Math.max(amin, squareMagnitude[i][j]));
                d -= 10 * Math.log10(Math.max(amin, squareMagnitude[i][j]));
                logSpec[i][j] = d;
                maxSpec = Math.max(maxSpec, d);
            }
        }

        for(int i=0;i<magnitude.size();i++) {
            for(int j=0;j>magnitude.get(0).size();j++) {
                logSpec[i][j] = Math.max(logSpec[i][j], maxSpec - topDb);
            }
        }

        xdb = logSpec;
    }

    /**
     * Calculate the fft and get magnitudes of signal
     * @param signal
     * @param numberofpoints
     * @return
     */
    public ArrayList<Complex> calculateFFT(float[] signal,int numberofpoints, int NFFT)
    {
        Complex[] y;
        Complex[] complexSignal = new Complex[numberofpoints];
        ArrayList<Double> absSignal = new ArrayList<Double>();

        //fill complex signal
        for(int i = 0; i < numberofpoints; i++){
            complexSignal[i] = new Complex(signal[i],0);
        }

        //do the fft
        y = new FFT().rfft(complexSignal, NFFT);

        ArrayList<Complex> ret = new ArrayList<>();

        assert  y.length == (NFFT / 2 + 1);

        for(int i=0;i<y.length;i++) {
            ret.add(y[i]);
        }

        //calculate magnitude and add to result
        // for(int i = 0; i < numberofpoints/2; i++){
        //     absSignal.add(i,Math.sqrt(Math.pow(y[i].re(), 2) + Math.pow(y[i].im(), 2)));
        // }
        return ret;
    }

}//end of class