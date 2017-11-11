// marr.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <iomanip>
#include <stdlib.h>
#include <sstream>
int inp[512][512];
int newinp[512][512];
int im[512][512];
using namespace std;

/*** New Data Types ***/
typedef unsigned char BYTE;

// create data structure to hold image
struct PIC
{
    unsigned int nChannel;
    bool InterLeaved;
    unsigned int Width, Height, Maxval;
    BYTE *img;
};

//function that loads in the header
bool LoadP5Header(ifstream &infile, PIC &pic)
{
    bool rtv = true;
    char buf[16];
    int bufIndex;
    int width, height, maxval;
    infile.read(buf, 2); // get the magic number
    buf[2]='\0';

    if(buf[0] == 'P' && buf[1] == '5')
	{
        infile.read(buf, 1);
        while(isspace(buf[0])) // Skip white space(s)
		{
            infile.read(buf,1);
		}

        // get width
        bufIndex = 0;
        while(bufIndex < 15 && !isspace(buf[bufIndex]))
		{
            bufIndex++;
            infile.read(buf+bufIndex, 1);
        }
        buf[bufIndex] = NULL;  // null terminated string
        width = atoi(buf);

        // get height
        infile.read(buf,1);
        while(isspace(buf[0])) // Skip white space(s)
		{
            infile.read(buf,1);
        }
        bufIndex = 0;
        while(bufIndex < 15 && !isspace(buf[bufIndex]))
		{
            bufIndex++;
            infile.read(buf+bufIndex, 1);
        }
        buf[bufIndex] = NULL;  // null terminated string
        height = atoi(buf);

       // get Maxval
		infile.read(buf,1);
        while(isspace(buf[0])) // Skip white space(s)
		{
            infile.read(buf,1);
        }
        bufIndex = 0;
        while(bufIndex < 15 && !isspace(buf[bufIndex]))
		{
            bufIndex++;
            infile.read(buf+bufIndex, 1);
        }
        buf[bufIndex] = NULL;  // null terminated string
        maxval = atoi(buf);

		// Skip white space(s)
		infile.read(buf,1);

        // set the image information in the struct
        pic.InterLeaved = false;
        pic.Width = width;
        pic.Height = height;
		pic.Maxval = maxval;

    }
    else rtv = false;

    return rtv;
}; // end of LoadP5Header()

//function that accepts an infile object and a PIC Object
//and reads in the PIC object
void LoadImage(ifstream &infile, PIC &pic)
{
    infile.read(reinterpret_cast<char *>(pic.img), pic.Width*pic.Height);
}

//function that accepts an outstream file and a PIC object to
//and writes the output stream to the PIC object
void WritePGM(ofstream & outfile, PIC pic)
{
    outfile << "P5" << endl;
    outfile << pic.Width << " " << pic.Height << endl;
	outfile << pic.Maxval << endl;

    outfile.write(reinterpret_cast<char *>(pic.img), pic.Width*pic.Height);
}


int _tmain(int argc, _TCHAR* argv[])
{
	
	
	ifstream infile_Pin;                /* input file */
    ofstream outfile_Pout;            /* output file */
    char inFileName_Pin[128]="lenna.pgm";            /* input file name */
    char outFileName_Pout[128]="marresult.pgm";        /* ouput file name */
    int scale;
	
	int Npixels, i,l;
	
    
    PIC Pin;            //source image
    PIC Pout;            //output file

    /* open input/output file */
    infile_Pin.open(inFileName_Pin, ios::binary);
    outfile_Pout.open(outFileName_Pout, ios::binary);

	/* check input file */
    if(!infile_Pin)
    {
        cout<<"Cannot open input file "<< inFileName_Pin <<endl;
        return 1;
    }

    if (LoadP5Header(infile_Pin, Pin)) // load pgm (Pin) image header information
    {
 
        Pin.img = new BYTE[Pin.Width*Pin.Height];

        LoadImage(infile_Pin, Pin);
	
	
        int w,h;
        w=(Pin.Width);
        h=(Pin.Height);
       
        	// allocate memory for output images
        Pout.nChannel = 1;
        Pout.InterLeaved = false;
        Pout.Width = w;
        Pout.Height = h;
		Pout.Maxval = Pin.Maxval;
		
		Npixels = w*h;//Pin.Width*Pin.Height;
        Pout.img = new BYTE [Npixels];
        int j=0,k;
    
       for (k=0;k<h;k++)
       {
          for (i=0;i<w;i++)
           {
            inp[k][i]=Pin.img[j];
			newinp[k][i]=0;
			im[k][i]=0;
			    j++;                 
           } 
       }

// LoG mask
	   int logmask[5][5]={{0, 0, 1, 0, 0},
	                  {0, 1, 2, 1, 0},
	                  {1, 2, -16, 2, 1},
	                  { 0, 1, 2, 1, 0},
	                  {0, 0, 1, 0, 0}};
	int l1,l2=3;
   int temp;l1=5;

j=0;l=0;


int border1=(l1/2);int border2=(l2/2);

// LoG convolution
int m;
int spix;
//
for (i=border1;i<(h-border1);i++)
{
    for (j=border1;j<(w-border1);j++)
	{
        spix=0;
        for ( k=0;k<l1;k++)
		{           
			for (m=0;m<l1;m++)
			{ 
				spix= spix+(logmask[k][m]*(inp[i-border1+k][j-border1+m]));}
		    }   
       
           newinp[i][j]=spix;
	}
}


 
// zero crossing & thresholding
int max=0, th;
for (k=0;k<h;k++)

{

for (i=0;i<w;i++)

{

if (newinp[k][i]>max)

     {max= newinp[k][i];}

}

} 
 
th=int(.075*max);

i=0;j=0;

//zero corossing

for (i=border2;i<(h-border2);i++)
{    j=0;
    for (j=border2;j<(w-border2);j++)
	{
        if ( newinp[i][j]!=0)
		{
			if ((newinp[i][j+1]>=0 && newinp[i][j-1]<0) || (newinp[i][j+1]<0 && newinp[i][j-1]>=0))
			      {    
					  if ((abs( newinp[i][j+1])-abs(newinp[i][j-1]))>0 && (newinp[i][j]>th))
			               { 
							   im[i][j]=255;
			                } 
			      }
            
			else if ((newinp[i+1][j]>=0 && newinp[i-1][j]<0) || (newinp[i+1][j]<0 && newinp[i-1][j]>=0))
			        { if ((abs(newinp[i+1][j])-abs(newinp[i-1][j]))>0 && newinp[i][j]>th)
			             {  im[i][j]=255;
			             }
			         }

			else if ((newinp[i+1][j+1]>=0 && newinp[i-1][j-1]<0) || (newinp[i+1][j+1]<0 && newinp[i-1][j-1]>=0))
			         { if ((abs(newinp[i+1][j+1])-abs(newinp[i-1][j-1]))>0 && newinp[i][j]>th)
			               { im[i][j]=255;
			               }
			         }
			
			else if ((newinp[i-1][j+1]>=0 && newinp[i+1][j-1]<0) || (newinp[i-1][j+1]<0 && newinp[i+1][j-1]>=0))
			        { if ((abs(newinp[i-1][j+1])-abs(newinp[i+1][j-1]))>0 && (newinp[i][j]>th))
			             { im[i][j]=255;
			             }
			        }
			
					}
	}
  }


j=(sizeof(inp)/2);


     
     m=0;
       for (k=0;k<h;k++)
       {
          for (i=0;i<w;i++)
           {
			   Pout.img[m]=im[k][i];
            m++;                       
           } 
       }            
          
     WritePGM(outfile_Pout, Pout);
		
	}

	return 0;
}

