//=================================================================================
//  Code for computing a new image quality metric
//
//
//  AUTORIGHTS
//  Copyright (C) 2017 Ecole Polytechnique Federale de Lausanne (EPFL), Switzerland.
//
//  Created by Radhakrishna Achanta on 07/02/2017.
//=================================================================================


#include<mex.h>
//#include <cstdio>
#include <cmath>
#include <cfloat>
#include <vector>
using namespace std;

struct RECT
{
    int x1, y1, x2, y2;
    double energy;
    int level;
    bool ignore;
};

//struct INTEGRAL
//{
//    INTEGRAL(const double* input, const int width, const int height)
//    {
//        ww = width+1;
//        hh = height+1;
//        int szz = ww*hh;
//        intimg.resize(szz,0);
//        for(int yy = 1, i = 0; yy < hh; yy++)
//        {
//            int ii = yy*ww;
//            for(int xx = 1; xx < ww; xx++, i++)
//            {
//                intimg[ii+xx] = intimg[ii+xx-1] + intimg[ii+xx-ww] - intimg[ii+xx-ww-1] + input[i];
//            }
//        }
//    }
//    double getSum(const int x1, const int y1, const int x2, const int y2) const
//    {
//        int x = x1; int xx = x2+1;
//        int y = y1; int yy = y2+1;
//        double val = intimg[yy*ww+xx] - intimg[y*ww+xx] - intimg[yy*ww+x] + intimg[y*ww+x];
//        return val;
//    }
//    
//    vector<double> intimg;
//    int ww, hh;
//};

//==============================================================================
///	maxMinDiff
//==============================================================================
double maxMinDiff(const RECT& r, double* energy, const int& width, const int& height)
{
    double maxval = 0;
    double minval = DBL_MAX;
    for(int y = r.y1; y < r.y2; y++)
    {
        for(int x = r.x1; x < r.x2; x++)
        {
            int i = y*width + x;
            if(maxval < energy[i]) maxval = energy[i];
            if(minval > energy[i]) minval = energy[i];
        }
    }
    return maxval-minval;
}

//==============================================================================
///	computeGradients
///
/// Assume 1-channel input.
//==============================================================================
void computeGradients(unsigned char* inp, const int width, const int height, int* edgemags)
{
    const int sz = width*height;
    const int w = width;
    
    for(int y = 1; y < height-1; y++)
    {
        for(int x = 1; x < width-1; x++)
        {
            int i = y*width+x;
            int l = i-1; int r = i+1;
            int t = i-w; int b = i+w;
            double dx = fabs(inp[l]-inp[r]);
            double dy = fabs(inp[t]-inp[b]);
            edgemags[i] = dx+dy;
        }
    }
}

//===============================================================================
// buildQuadtree
//===============================================================================
void buildQuadtree(double* energy, const int width, const int height, vector<RECT>& rects)
{
    RECT tr;
    tr.x1 = 0;tr.x2 = width;tr.y1 = 0;tr.y2 = height;
    tr.level = 0;
    tr.energy = maxMinDiff(tr,energy,width,height);//intimg.getSum(tr.x1, tr.y1, tr.x2-1, tr.y2-1);
    tr.ignore = true;
    rects.push_back(tr);
    int rectcount = 1;
    for( int j = 0; j < rectcount ; j++ )
    {
        if( min( rects[j].x2-rects[j].x1, rects[j].y2-rects[j].y1) > 4 )//if either side is > 4 split.
        {
            bool condition = rects[j].energy > 1;
            if( true == condition )
            {
                rects[j].ignore = true;
                //getFourQuadrants(rects, j, rects[j].level+1, intimg);
                //---------------------------------------------------------------------
                vector<RECT> rect(4);
                RECT r = rects[j];
                double halfx = (r.x1+r.x2)/2;
                double halfy = (r.y1+r.y2)/2;
                
                int ri = 0;
                rect[ri+0].x1 = r.x1;
                rect[ri+0].y1 = r.y1;
                rect[ri+0].x2 = halfx;
                rect[ri+0].y2 = halfy;
                
                rect[ri+1].x1 = halfx;
                rect[ri+1].y1 = r.y1;
                rect[ri+1].x2 = r.x2;
                rect[ri+1].y2 = halfy;
                
                rect[ri+2].x1 = r.x1;
                rect[ri+2].y1 = halfy;
                rect[ri+2].x2 = halfx;
                rect[ri+2].y2 = r.y2;
                
                rect[ri+3].x1 = halfx;
                rect[ri+3].y1 = halfy;
                rect[ri+3].x2 = r.x2;
                rect[ri+3].y2 = r.y2;
                
                for(int n = 0; n < 4; n++)
                {
                    rect[n].ignore = false;
                    rect[n].level = rects[j].level+1;
                    rect[n].energy = maxMinDiff(rect[n],energy,width,height);
                    rects.push_back(rect[n]);
                }
                //---------------------------------------------------------------------
                rectcount += 4;
            }
        }
    }
}

//===============================================================================
// Inpur argument order must be reference image followed by candidate image
//
// [7 FEb 2017]
//===============================================================================
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 2)
    {
        mexErrMsgTxt("The correct number of arguments is two, namely, the refence image and the candidate image") ;
    }
    //---------------------------
    const int numelements   = mxGetNumberOfElements(prhs[0]) ;
    const mwSize numdims = mxGetNumberOfDimensions(prhs[0]);
    const mwSize* dims  = mxGetDimensions(prhs[0]) ;
    unsigned char* refbytes = (unsigned char*)mxGetData(prhs[0]) ;//mxGetData returns a void pointer, so cast it
    unsigned char* imgbytes = (unsigned char*)mxGetData(prhs[1]) ;
    const int width = dims[1];
    const int height = dims[0];//Note: first dimension provided is height and second is width
    const int sz = width*height;
    const int szsz = sz+sz;
    
    //---------------------------
    // Allocate memory
    //---------------------------
    unsigned char* rref  = (unsigned char*)mxMalloc( sizeof(unsigned char)      * sz );
    unsigned char* gref  = (unsigned char*)mxMalloc( sizeof(unsigned char)      * sz );
    unsigned char* bref  = (unsigned char*)mxMalloc( sizeof(unsigned char)      * sz );
    unsigned char* rinp  = (unsigned char*)mxMalloc( sizeof(unsigned char)      * sz );
    unsigned char* ginp  = (unsigned char*)mxMalloc( sizeof(unsigned char)      * sz );
    unsigned char* binp  = (unsigned char*)mxMalloc( sizeof(unsigned char)      * sz );
    
    int* reng  = (int*)mxMalloc( sizeof(int)      * sz );
    int* geng  = (int*)mxMalloc( sizeof(int)      * sz );
    int* beng  = (int*)mxMalloc( sizeof(int)      * sz );
    double* energy_ref  = (double*)mxMalloc( sizeof(double)      * sz );
    
    for(int i = 0; i < sz; i++)
    {
        reng[i] = geng[i] = beng[i] = 0;
    }

    
    //---------------------------
    // Read the pixel values
    //---------------------------
    //if(2 == numdims)
    if(numelements/sz == 1)//if it is a grayscale image, copy the values directly into the lab vectors
    {
        for(int x = 0, ii = 0; x < width; x++)//reading data from column-major MATLAB matrics to row-major C matrices (i.e perform transpose)
        {
            int ind = 0;
            for(int y = 0; y < height; y++)
            {
                //int i = y*width+x;
                int i = ind + x;
                rref[i] = refbytes[ii];
                gref[i] = refbytes[ii];
                bref[i] = refbytes[ii];
                rinp[i] = imgbytes[ii];
                ginp[i] = imgbytes[ii];
                binp[i] = imgbytes[ii];
                ii++;
                ind += width;
            }
        }
    }
    else//else convert from rgb to lab
    {
        for(int x = 0, ii = 0; x < width; x++)//reading data from column-major MATLAB matrics to row-major C matrices (i.e perform transpose)
        {
            int ind = 0;
            for(int y = 0; y < height; y++)
            {
                //int i = y*width+x;
                int i = ind + x;
                rref[i] = refbytes[ii];
                gref[i] = refbytes[ii+sz];
                bref[i] = refbytes[ii+szsz];
                rinp[i] = imgbytes[ii];
                ginp[i] = imgbytes[ii+sz];
                binp[i] = imgbytes[ii+szsz];
                ii++;
                ind += width;
            }
        }
    }

    
    computeGradients(rref,width,height,reng);
    computeGradients(gref,width,height,geng);
    computeGradients(bref,width,height,beng);
    for(int i = 0; i < sz; i++)
    {
         energy_ref[i] = (reng[i]+geng[i]+beng[i])/3.0;
    }
    vector<RECT> qtree(0);
    buildQuadtree(&energy_ref[0],width,height,qtree);
    
    //-----------------------------------
    // Compute mean error over the quadtree
    //-----------------------------------
    double meanerrsum = 0;
    if(1)//quadtree based mse
    {
        int validrectcount = 0;
        int numrects = qtree.size();
        for( int n = 0; n < numrects; n++ )
        {
            if( false == qtree[n].ignore )
            {
                validrectcount++;
                int x1 = qtree[n].x1; int y1 = qtree[n].y1; int x2 = qtree[n].x2; int y2 = qtree[n].y2;
                double area = (x2-x1)*(y2-y1);
                double rrsum(0),risum(0),grsum(0),gisum(0),brsum(0),bisum(0);
                for( int y = qtree[n].y1; y < qtree[n].y2; y++ )
                {
                    for( int x = qtree[n].x1; x < qtree[n].x2; x++ )
                    {
                        int i = y*width+x;
                        rrsum += rref[i]; risum += rinp[i];
                        grsum += gref[i]; gisum += ginp[i];
                        brsum += bref[i]; bisum += binp[i];
                    }
                }
                double rerr = (rrsum-risum)/area;
                double gerr = (grsum-gisum)/area;
                double berr = (brsum-bisum)/area;
                meanerrsum += rerr*rerr+gerr*gerr+berr*berr;
            }
        }
        meanerrsum /= (3.0*validrectcount);
    }
    double qpsnr = 10.0*log10((255.0*255.0)/meanerrsum);
    
//    FILE* pf = fopen("/Users/radhakrishnaachanta/rktemp/qpsnr.txt","w");
//    fprintf(pf,"qpsnr = %f\n", qpsnr);
//    fclose(pf);
    
    //---------------------------
    // Assign the QPSNR value
    //---------------------------
    plhs[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
    double* qpsnrval = (double*)mxGetData(plhs[0]);//gives a void*, cast it to int*
    *qpsnrval = qpsnr;
    
    //---------------------------
    // Deallocate memory
    //---------------------------
    mxFree(rref);
    mxFree(gref);
    mxFree(bref);
    mxFree(rinp);
    mxFree(ginp);
    mxFree(binp);
    mxFree(reng);
    mxFree(geng);
    mxFree(beng);
    mxFree(energy_ref);
}


