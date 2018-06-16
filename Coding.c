#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


#define height 512
#define width 512
#define m 2
#define max_times 300
#define delta 0.00001
#define cluster 128
#define MAX 256
#define e 2.718281828459
#define pi 3.1415926535897932384626

/************************ Global Data ******************************/
double write[height][width], write1[height][width], write2[height][width];
double J[max_times];
double membership[cluster][MAX];
double distance[cluster][MAX];
double center[cluster];
int freq[MAX];
int Function_check=0;

/*************** The Function of Fuzzy & mapping *******************/
void Initialize_Membership_Degree();
void Fuzzy_C_Mean();
void get_ClusterCenter(int t);
void Update_Membership_Degree(int t);
void find_Cluster();

/****************** objective quantitative measurement **********************/
void ABME();
void ENTROPY();
void EME();
void MSE();
void PSNR();

void data_read()
{
	int i,j;
	FILE *org;
    org =  fopen("C://Users//Administrator//Desktop//FuzzyC_LinearScaleUp//ct_2.raw", "rb");
    /**  fopen("c://RAW_file//MRI_h.raw", "rb"); **/

   for (i = 0; i < height; i++)
   {
      for (j = 0; j < width; j++)
      {
         write[i][j] = fgetc(org);
      }
   }
   printf("Original Image(8x8): \n");
   for (i = 0; i < 8; i++)
   {
      for (j = 0; j < 8; j++)
      {
         printf("%f\t",write[i][j]);
      }
	   printf("\n");
   }
   fclose(org);



}

void ABME()
{
	int i,j,k;
	double write_avg=0;
	double write1_avg=0;
	double abme=0;

	for(i=0; i<height; i++)
	{
		for(j=0; j<width; j++)
		{
			write_avg += write[i][j];
			write1_avg += write1[i][j];
		}
	}
	write_avg = write_avg/(height*width);
	write1_avg = write1_avg/(height*width);
	abme = fabs(write_avg - write1_avg);
	printf("ABME = %f\n", abme);
}

void ENTROPY(double img[height][width])
{
	int i,j,k;
	int freq1[MAX]={0};
	double pr;
	double entropy=0;

	for(k=0; k<MAX; k++)
	{
		for(i=0; i<height; i++)
		{
			for(j=0; j<width; j++)
			{
				if(img[i][j]==k)
				{
					freq1[k]++;
				}
			}
		}
	}
	for(k=0; k<MAX; k++)
	{
		pr = (double)freq1[k]/(height*width);
		if(pr != 0) entropy += -pr*log(pr)/log(2.0);
	}
	printf("ENTROPY = %f\n", entropy);
}

void EME(double img[height][width])
{
	int h,w,i,j;
	int L=5;
	int how_many = floor((double)height/L);
	int block_min, block_max;
	double block_ratio=0;
	double eme=0;

	for(h=0; h<how_many; h++)
	{
		for(w=0; w<how_many; w++)
		{
			block_min=1000; block_max=-1;
			for(i=h*L; i<h*L+5; i++)
			{
				for(j=w*L; j<w*L+5; j++)
				{
					if(img[i][j]<=block_min)
						block_min = img[i][j];
					if(img[i][j]>=block_max)
						block_max = img[i][j];
				}
			}
			if(block_min>0)
			{
				block_ratio = (double)block_max/block_min;
				eme += 20.0*log(block_ratio)/log(10.0);
			}
		}
	}
	eme = eme/(how_many*how_many);
	printf("EME = %f\n", eme);
}

void MSE()
{
	int i,j;
	double mse=0;

	for(i=0; i<height; i++)
    {
        for(j=0; j<width; j++)
        {
            mse += pow((write1[i][j]-write[i][j]), 2)/(height*width);
        }
    }
	printf("MSE = %f\n", mse);
	PSNR(mse);
}

void PSNR(double mse)
{
	int i,j;
	double psnr=0;

	psnr = 10.0*log((MAX-1)*(MAX-1)/mse)/log(10.0);
	printf("PSNR = %f\n", psnr);
}


int main()
{
    int i,j, k;
	clock_t start, finish; //set time
	double totaltime; //set time

	FILE *fuzzy_c_means;
    fuzzy_c_means = fopen("d://Medical Image//resImg.raw", "wb");
	start=clock(); //set time
	
    data_read();
	/***************************** Image Histogram *********************************/
    for(k=0; k<MAX; k++)
    {      freq[k] = 0;
        for(i=0; i<height; i++)
        {	for(j=0; j<width; j++)
            {	
				if(write[i][j]==k) 
					 freq[k]++;	
		    }
        }
    }
    //for(k=0; k<MAX; k++){	printf("freq[%d]=%d\n", k, freq[k]);	}


    /********************** Image Fuzzy C-Means cluster process ***************************/
    Initialize_Membership_Degree();/** set chushihua membership degree **/
    Fuzzy_C_Mean();/** implement Fuzzy C-Means Function **/

	/**************************** Image Data Output ***********************************/
	for(i=0; i<height; i++)
    {
        for(j=0; j<width; j++)
        {
            fputc(write1[i][j], fuzzy_c_means);
        }
    }
	/************* Test ABME **************/
     ABME();
	/************ Test ENTROPY ************/
	 ENTROPY(write); ENTROPY(write1);
	/************* Test EME ***************/
	 EME(write);  EME(write1);
	/********** Test MSE & PSNR ***********/
	 MSE();

	finish=clock(); //set time
    totaltime=(double)(finish-start)/CLOCKS_PER_SEC; //set time
    printf("\nthis program running time: %f s\n", totaltime); //set time
    return 0;
}

/** set chushihua membership degree and make their sum is one **/
void Initialize_Membership_Degree()
{
    int i, j, k;
    double member_sum;

    for(j=0; j<MAX; j++)
    {
        member_sum = 0;
        for(i=0; i<cluster; i++)
        {
            membership[i][j] = rand();/** 0~32767 rand numbers **/
            member_sum = member_sum + membership[i][j];
        }
        for(i=0; i<cluster; i++)
        {
            membership[i][j] = membership[i][j]/member_sum;/**  Membership_degree of their Sum is One **/
        }

    }

}

void Fuzzy_C_Mean()
{
    int i, t;

    for(t=1; t<=max_times; t++)
    {
        get_ClusterCenter(t);/** get one Center of Image Cluster, all 4 cluster center **/
        Update_Membership_Degree(t);/** get new membership degree **/
        if(t>1)
            if(fabs(J[t]-J[t-1]) < delta || Function_check==1)
                break;
    }

    printf("\n");
    
    find_Cluster();/** find cluster of every pixel through the Max_membership **/


}

/** calculate center of cluster **/
void get_ClusterCenter(int t)
{
    int i, j;
    double fenZ, fenM;

    for(i=0; i<cluster; i++)
    {
        fenZ=0; fenM=0;
        for(j=0; j<MAX; j++)
        {
            fenZ += pow(membership[i][j], m)*j*freq[j];
            fenM += pow(membership[i][j], m)*freq[j];
        }
        center[i] = fenZ/fenM;
    }

}

/** get new membership degree matrix **/
void Update_Membership_Degree(int t)
{
    int i, j, k, n;
    double sum_dis;
    double member_sum;
    J[t]=0;

    /** get Euclid equation and Object function **/
    for(j=0; j<MAX; j++)
    {
        for(i=0; i<cluster; i++)
        {
			distance[i][j] = fabs(j-center[i]); 
			if(distance[i][j]==0){ 
				Function_check=1; 
				return;
			}
            J[t] +=  pow(membership[i][j], m)*pow(distance[i][j], 2)*freq[j];
        }
    }

    /** get New membership_degree **/
    for(j=0; j<MAX; j++)
    {
        for(i=0; i<cluster; i++)
        {
            sum_dis = 0;
            for(k=0; k<cluster; k++)
            {
                sum_dis += pow(distance[k][j], -2/(m-1));//sum_dis += pow(distance[i][j]/distance[k][j], 2/(m-1));
            }
            membership[i][j] = pow(distance[i][j], -2/(m-1))/sum_dis; //membership[i][j] = 1/sum_dis;
        }
    }

    /** get membership_degree_sum to be One **/
    for(j=0; j<MAX; j++)
    {
        member_sum = 0;
        for(i=0; i<cluster; i++)
        {
            member_sum += membership[i][j];
        }
        for(i=0; i<cluster; i++)
        {
            membership[i][j] = membership[i][j]/member_sum;/**  Membership_degree of their Sum is One **/
        }
    }

    printf("t=%d\tJ[%d]=%f\n", t, t, J[t]);

}

/** find cluster that every pixel belong to through the Max_membership(possibility is big) **/
void find_Cluster()
{

    int i, j, k, t, time;
    int classNum;
    double max_center;  double min_center;
    double max_memDegree;
    int tf_center[cluster]={0};
	double order_C128[cluster]={0};
    double order_C256[MAX]={0};
    int index_C256[MAX]={0};
    int compare_center[cluster]={0};
    int target_center[MAX]={0};
	double mean_C=0, sigma_C=0;
	double R[cluster];// results by Linear scaling-up function

	//************** get mean and sigma based on the 128 center values *************//
	for(i=0; i<cluster; i++)
	{
		mean_C += center[i];
	}
	mean_C = mean_C/cluster;

	for(i=0; i<cluster; i++)
	{
		sigma_C += fabs(center[i] - mean_C);
	}
	sigma_C = sigma_C/cluster;
	//printf("mean_C = %f, sigma_C = %f\n", mean_C, sigma_C);

	//*************** R[i] by Linear scaling-up function ***********************//
	for(i=0; i<cluster; i++)
	{
		R[i] = (MAX-1) * pow((double)i/cluster, mean_C/sigma_C); //printf("R[%d]=%f\n", i,R[i]);
	}

    //*********** get 128 cluster center in the order of small to large ***********//
   for(i=0; i<cluster; i++)
    {
        max_center = 1000000;
        for(j=0; j<cluster; j++)
        {
            if(tf_center[j]==0 && center[j]<max_center)
            {
                max_center = center[j];
                k=j;
            }
        }
        tf_center[k]=1; order_C128[i] = max_center;
    }
    //printf("\n"); for(i=0; i<cluster; i++){ printf("order_C128[%d] = %f\n", i, order_C128[i]); } printf("\n\n");

    //*************** compare order_C128[i] and R[i] ***********************//
    for(i=0; i<cluster; i++)
    {
        if(order_C128[i] > R[i])
            compare_center[i] = (int)(order_C128[i] + 0.5);
        else
            compare_center[i] = (int)(R[i] + 0.5);
    }


    //************* jth gray value -> center value of ith cluster ****************//
    for(j=0; j<MAX; j++)
    {
        max_memDegree=-1;
        for(i=0; i<cluster; i++)
        {
            if(membership[i][j]>max_memDegree)
            {
                max_memDegree = membership[i][j];
                classNum = i;
            }
        }
        order_C256[k] = center[classNum];
    }

    //************* connect order_C128 and order_C256 *************************//
    for(i=0; i<cluster; i++)
    {
        for(j=0; j<MAX; j++)
        {
            if(order_C256[j] == order_C128[i])
                index_C256[j] = i;
            else
                continue;
        }
    }
    for(i=0; i<cluster; i++)
    {
        for(j=0; j<MAX; j++)
        {
            if(index_C256[j] == i)
                target_center[j] = compare_center[i];
        }
    }

    //********* obtained target_center send to write1[i][j] matrix **********//
    for(k=0; k<MAX; k++)
    {
        for(i=0; i<height; i++)
        {
            for(j=0; j<width; j++)
            {
                if(write[i][j]==k)
                {
                    write1[i][j] = target_center[k];
                }

            }
        }
    }

}
