/*
Open-Source AutoShim Copyright notice
================
 
The AutoShim source is copyrighted, but you can freely use and copy it
as long as you don't change or remove the copyright notice:
 
----------------------------------------------------------------------
Copyright (c) 2007, Novartis Vaccines and Diagnostics, Inc., a Delaware corporation with its principal offices at 4560 Horton Street, Emeryville, CA 94530

All rights reserved.
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
·	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer. 
·	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution. 
·	Neither the name of the Novartis Vaccines and Diagnostics, Inc. nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission. 
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Written by David C. Sullivan, 2007

compile by:
gcc -c RemoveDuplicatesFromSD_9.c ; gcc -lm -g RemoveDuplicatesFromSD_9.o -o RemoveDuplicatesFromSD_9

usage:

*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXSQUARE 1
#define MAXLINESRECORD 1000
#define MAXCONFS 3000

int FirstPassFilter(char ***confs,int confcnt,int confnumlines);

int main(int argc, char *argv[])
{
char line[10240];
char CurrentChID[128],id[128],dollars[128];
char ***CurrentMol; /* MAXCONFS confs, MAXLINESRECORD lines each, 100 characters wide */
int end, linecount,linecountFirst, confcount=0,i,j,firsttimeflag=1;
FILE *infile;

	/* allocate memory for conformations of a chironid */
	CurrentMol = (char ***)calloc(MAXCONFS,sizeof(char **));
	for(i=0;i<MAXCONFS;i++)
	{
		CurrentMol[i] = (char **)calloc(MAXLINESRECORD,sizeof(char *));
		for(j=0;j<MAXLINESRECORD;j++)
			CurrentMol[i][j] = (char *)calloc(100,sizeof(char));
	}
	

	infile = fopen(argv[1], "r");
	while (fgets(line, sizeof(line), infile) != NULL)
	{

		sscanf(line,"%s",id);
   /*             print("%s",id);  */
		if(strcmp(id,CurrentChID) != 0)
		{
                

			if(confcount > 1)
			{

				FirstPassFilter(CurrentMol,confcount,linecountFirst);
			}


			strcpy(CurrentChID,id);
                    /*    print("%s",CurrentChID); */
	            		if(strlen(CurrentChID) != 27)
				fprintf(stderr,"Error: Current Chironid of wrong length: %d   %s\n",strlen(CurrentChID),CurrentChID);

			strcpy(CurrentMol[1][1],line);

			confcount = 1;
			end = 0;
			linecount = 1;
			while(end == 0)
			{
				linecount++;
				fgets(line, sizeof(line), infile);
				strcpy(CurrentMol[confcount][linecount],line);
				sscanf(line,"%s",dollars);
				if(strcmp(dollars,"$$$$") == 0)
					end = 1;
			}
			fgets(line, sizeof(line), infile);
			linecountFirst = linecount;
		}
		else
		{

			confcount++;
			end = 0;

			if(confcount>2995)fprintf(stderr,"Error: too many conformations: %d  %s\n",confcount,id);			

			strcpy(CurrentMol[confcount][1],line);

			linecount = 1;
			while(end == 0)
			{
				linecount++;
				fgets(line, sizeof(line), infile);
				strcpy(CurrentMol[confcount][linecount],line);
				sscanf(line,"%s",dollars);
				if(strcmp(dollars,"$$$$") == 0)
					end = 1;
			}
			fgets(line, sizeof(line), infile);
			if(linecount != linecountFirst)
				fprintf(stderr,"Error: linecounts within chironid not equal %s %s   %d %d\n",id,CurrentChID,linecount,linecountFirst);
		}
	}
	fclose(infile);


	FirstPassFilter(CurrentMol,confcount,linecountFirst);


	
return 0;
}



int FirstPassFilter(char ***confs,int confcnt,int confnumlines)
{
int i,j,m,n,d,cnt,heavyatomcount,*keeplist;
double ***confxyz,maxdiffSq,diffSq,**confdistmatrix,drmsd;
double maxX,maxY,maxZ;

	confxyz = (double ***)calloc(confcnt+1,sizeof(double **));
	keeplist = (int *)calloc(confcnt+1,sizeof(int));
	for(m=0;m<=confcnt;m++)
	{
		keeplist[m] = 1;
		confxyz[m] = (double **)calloc(200,sizeof(double *));
		for(i=0;i<200;i++)
			confxyz[m][i] = (double *)calloc(4,sizeof(double));
	}
	
	
	for(m=1;m<=confcnt;m++)
	{
		heavyatomcount = 0;
		for(i=0;i<confnumlines;i++)
		{
			if((confs[m][i][31] != 'H')&&(confs[m][i][35] != ' ') && (strlen(confs[m][i]) >= 42))
			{
				heavyatomcount++;
				sscanf(confs[m][i],"%lf %lf %lf",
				&confxyz[m][heavyatomcount][1],&confxyz[m][heavyatomcount][2],&confxyz[m][heavyatomcount][3]);
			}
		}
	}





	for(m=1;m<=confcnt;m++)
	{
		if(keeplist[m] == 1)
		{
		for(n=m+1;n<=confcnt;n++)
		{
		
			/* max difference measure */
			maxdiffSq = 0;

			for(i=1;i<=heavyatomcount;i++)
			{
			     diffSq =	((confxyz[m][i][1] - confxyz[n][i][1]) * (confxyz[m][i][1] - confxyz[n][i][1])) + 
					((confxyz[m][i][2] - confxyz[n][i][2]) * (confxyz[m][i][2] - confxyz[n][i][2])) + 
					((confxyz[m][i][3] - confxyz[n][i][3]) * (confxyz[m][i][3] - confxyz[n][i][3]));

				if(diffSq > maxdiffSq)
				{
					maxdiffSq = diffSq;

					maxX = confxyz[m][i][1];
					maxY = confxyz[m][i][2];
					maxZ = confxyz[m][i][3];
				}
			}
			if (maxdiffSq < MAXSQUARE)
			{
				keeplist[n] = 0;
/*keeplist[n] = 0 means that the pose is disgarded.*/
			}

		}
		}
	}
	
	cnt=0;
	for(m=1;m<=confcnt;m++)
	{	
		if(keeplist[m] == 1)
		{
			cnt++;
/*cnt is the count of the number of poses that have been processed. cnt<10 is the current setting meaning */
/*that a maximum of 9 conformations will be processed.*/
			if(cnt<16) 
			{
				for(i=1;i<=confnumlines;i++)
				{
					printf("%s",confs[m][i]);

				}
			}
		}
	}

return 0;
}


			
			
