package cn.hzw.doodle.STFT;

import android.util.Log;

public class FFT {

    /**
     * Compute the FFT of x[], assuming its length is a power of 2
     * @param x
     * @return
     */
    public Complex[] fft(Complex[] x) {
        int N = x.length;

        // base case
        if (N == 1) return new Complex[] { x[0] };

        // radix 2 Cooley-Tukey FFT
        if (N%2 != 0) {
            throw new RuntimeException("N is not a power of 2");
        }

        // fft of even terms
        Complex[] even = new Complex[N/2];
        for (int k = 0; k < N/2; k++) {
            even[k] = x[2*k];
        }
        Complex[] q = fft(even);

        // fft of odd terms
        Complex[] odd  = even;  // reuse the array
        for (int k = 0; k < N/2; k++) {
            odd[k] = x[2*k + 1];
        }
        Complex[] r = fft(odd);

        // combine
        Complex[] y = new Complex[N];
        for (int k = 0; k < N/2; k++) {
            double kth = -2 * k * Math.PI / N;
            Complex wk = new Complex(Math.cos(kth), Math.sin(kth));
            y[k]       = q[k].plus(wk.times(r[k]));
            y[k + N/2] = q[k].minus(wk.times(r[k]));
        }
        return y;
    }

    public Complex[] rfft(Complex[] x, int NFFT) {
        if(NFFT != 254) {
            throw new RuntimeException("nfft is not 254");          // 针对于NFFT是254的情况
        }
        int LENGTH = 256;
        Complex[] trunc = new Complex[LENGTH];

        for(int i=0;i<NFFT;i++) {
            trunc[i+1] = new Complex(x[i].re(), x[i].im());
        }

        trunc[0] = new Complex(trunc[1].re(), trunc[1].im());
        trunc[255] = new Complex(trunc[254].re(), trunc[254].im());

        Complex[] res = this.fft(trunc);
        Complex[] ret = new Complex[NFFT/2+1];

        assert res.length == NFFT;

        for(int i=0;i<NFFT/2+1;i++) {
            ret[i] = res[i+1];
        }

        return ret;
    }

    /**
     * Compute the inverse FFT of x[], assuming its length is a power of 2
     * @param x
     * @return
     */
    public Complex[] ifft(Complex[] x) {
        int N = x.length;
        Complex[] y = new Complex[N];

        // take conjugate
        for (int i = 0; i < N; i++) {
            y[i] = x[i].conjugate();
        }

        // compute forward FFT
        y = fft(y);

        // take conjugate again
        for (int i = 0; i < N; i++) {
            y[i] = y[i].conjugate();
        }

        // divide by N
        for (int i = 0; i < N; i++) {
            y[i] = y[i].times(1.0 / N);
        }

        return y;

    }

    /**
     * Compute the circular convolution of x and y
     * @param x
     * @param y
     * @return
     */
    public Complex[] cconvolve(Complex[] x, Complex[] y) {

        // should probably pad x and y with 0s so that they have same length
        // and are powers of 2
        if (x.length != y.length) { throw new RuntimeException("Dimensions don't agree"); }

        int N = x.length;

        // compute FFT of each sequence
        Complex[] a = fft(x);
        Complex[] b = fft(y);

        // point-wise multiply
        Complex[] c = new Complex[N];
        for (int i = 0; i < N; i++) {
            c[i] = a[i].times(b[i]);
        }

        // compute inverse FFT
        return ifft(c);
    }


    /**
     * Compute the linear convolution of x and y
     * @param x
     * @param y
     * @return
     */
    public Complex[] convolve(Complex[] x, Complex[] y) {
        Complex ZERO = new Complex(0, 0);

        Complex[] a = new Complex[2*x.length];
        for (int i = 0;        i <   x.length; i++) a[i] = x[i];
        for (int i = x.length; i < 2*x.length; i++) a[i] = ZERO;

        Complex[] b = new Complex[2*y.length];
        for (int i = 0;        i <   y.length; i++) b[i] = y[i];
        for (int i = y.length; i < 2*y.length; i++) b[i] = ZERO;

        return cconvolve(a, b);
    }

    /**
     * Display an array of Complex numbers to standard output
     * @param x
     * @param title
     */
    public void show(Complex[] x, String title) {
        System.out.println(title);
        System.out.println("-------------------");
        for (int i = 0; i < x.length; i++) {
            System.out.println(x[i]);
        }
        System.out.println();
    }

}//end of class

