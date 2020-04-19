package cn.hzw.doodle;

import android.content.Context;
import android.graphics.Canvas;
import android.graphics.Color;
import android.graphics.Paint;
import android.graphics.Path;
import android.graphics.Picture;
import android.graphics.Point;
import android.util.AttributeSet;
import android.view.View;

public class MathView extends View {
    private Point mCoo = new Point(500,700);
    private Picture mCooPicture;
    private Picture mGrindPicture;
    private Paint mHelpPrint;
    private Paint mPaint;
    private Path mPath;

    public MathView(Context context) {
        this(context,null);
    }

    public MathView(Context context, AttributeSet attrs) {
        super(context, attrs);
        init();
    }

    private void init() {
        mPaint = new Paint(Paint.ANTI_ALIAS_FLAG);
        mPaint.setColor(Color.BLUE);
        mPaint.setStrokeWidth(2);
        mPaint.setStyle(Paint.Style.STROKE);
        mPaint.setStrokeCap(Paint.Cap.ROUND);

        mPath = new Path();
    }

    @Override
    protected void onDraw(Canvas canvas) {

        for(int i=0;i<mydata.length-1;i++) {
            canvas.drawLine(i+180,mydata[i]/2+250,i+181,mydata[i+1]/2+250,mPaint);
        }

    }

    float[] mydata = new float[640];

    public void setData(float[] data) {
        mydata = data;
        invalidate();
    }

}

