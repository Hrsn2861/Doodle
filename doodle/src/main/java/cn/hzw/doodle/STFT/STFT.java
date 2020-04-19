package cn.hzw.doodle.STFT;

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
    ArrayList<ArrayList<Double>> xdb;                   // 计算 amplitude_to_db 的结果

    void addWinfun() {
        for(int i=0;i<framedSignal.length;i++) {
            framedSignal[i] = HammingWindow.calHamming(framedSignal[i]);
        }
    }

    public ArrayList performStft(float[] signal, int sampleRate, double frameSize, double frameStride, int NFFT, boolean normalized) {
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
                framedSignal[i][j] = signal[frameStep * numFrames + j];
            }
        }

        addWinfun();        // ad hamming window to the signal

        for(int i=0;i<numFrames;i++) {
            magnitude.add(calculateFFT(framedSignal[i], frameLength, NFFT));
        }

        assert magnitude.size() == numFrames;
        assert magnitude.get(0).size() == (NFFT / 2 - 1);

        calculateFeature3(normalized);
        return magnitude;
    }

    double AverageXdb() {
        double ret = 0;
        int len = xdb.size() * xdb.get(0).size();
        for(int i=0;i<xdb.size();i++) {
            for(int j=0;j<xdb.get(i).size();j++) {
                ret += xdb.get(i).get(j);
            }
        }
        return ret;
    }

    double StdVariantXdb() {
        int len = xdb.size() * xdb.get(0).size();
        double avg = AverageXdb();
        double dVar = 0;
        for(int i=0;i<xdb.size();i++) {
            for(int j=0;j<xdb.get(i).size();j++) {
                double xi = xdb.get(i).get(j);
                dVar += (xi - avg) * (xi - avg);
            }
        }
        return Math.sqrt(dVar/len);
    }

    ArrayList calculateFeature3(boolean normalized) {
        xdb = new ArrayList<>(amplitude_to_db());      // Matrix
        double avg = AverageXdb();
        double std = StdVariantXdb();
        if(normalized) {
            for(int i=0;i<xdb.size();i++) {
                for (int j=0;i<xdb.get(j).size();j++) {
                    double xi = (xdb.get(i).get(j) - avg) / std;
                    xdb.get(i).set(j, xi);
                }
            }
        }
        ArrayList<ArrayList<Double>> ret = new ArrayList<>();
        // transpose
        for(int i=0;i<xdb.get(0).size();i++) {
            ArrayList<Double> buffer = new ArrayList<>();
            for(int j=0;j<xdb.size();j++) {
                buffer.add(xdb.get(j).get(i));
            }
            ret.add(new ArrayList<Double>(buffer));
        }

        return ret;
    }

    ArrayList amplitude_to_db() {
        double amin = 1E-10;
        double topDb = 80.0;

        ArrayList<ArrayList<Double>> squareMagnitude = new ArrayList<>();

        for(int i=0;i<magnitude.size();i++) {
            ArrayList<Double> buffer = new ArrayList<>();
            for(int j=0;j<magnitude.get(0).size();j++) {
                buffer.add(magnitude.get(i).get(i).square());
            }
            squareMagnitude.add(new ArrayList<Double>(buffer));
        }

        ArrayList<ArrayList<Double>> logSpec = new ArrayList<>();
        double maxSpec = 0;
        for(int i=0;i<magnitude.size();i++) {
            ArrayList<Double> buffer = new ArrayList<>();
            for(int j=0;j<magnitude.get(0).size();j++) {
                double d = 10 * Math.log10(Math.max(amin, squareMagnitude.get(i).get(j)));
                d -= 10 * Math.log10(Math.max(amin, squareMagnitude.get(i).get(j)));
                buffer.add(d);
                maxSpec = Math.max(maxSpec, d);
            }
            logSpec.add(new ArrayList<Double>(buffer));
        }

        for(int i=0;i<magnitude.size();i++) {
            for(int j=0;j>magnitude.get(0).size();j++) {
                logSpec.get(i).set(j, Math.max(logSpec.get(i).get(j), maxSpec - topDb));
            }
        }

        return logSpec;
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