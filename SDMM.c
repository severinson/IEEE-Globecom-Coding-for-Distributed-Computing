#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

int main()
{
    int n,m,P,Pc,delta,t,s,d,ts,ss,ds;

    m = 36;
    n = 36;
    P = 3000;
    Pc = 0;
    
    double *threshold = new double[MIN(m,n)+1];
    double *load = new double[MIN(m,n)+1];
    
    memset(threshold,0,sizeof(double)*(MIN(m,n)+1));
    memset(load,0,sizeof(double)*(MIN(m,n)+1));

    int counter = 0;
    
    if (Pc >= P)
        return 0;
    
    for (int s=1;s<=MIN(n,m);s++)
        if ((m%s)==0 && (n%s)==0)

        {
            t = m / s;
            d = n / s;
            
            if (s < t)
            {
                if ((Pc %s)==0)
                    delta = Pc / s;
                else
                    delta = ((int) (Pc / s))+1;
                
                ts = t + delta;
                ds = d + delta;
                
                if (Pc==0)
                    threshold[counter++] = t*s*d + s -1;
                else if (Pc >= 1 && (Pc%s==0) && delta == Pc/s )
                    threshold[counter++] = ts*s*(d+1) + s*delta-1;
                else if (Pc >= 1 && (double) delta > Pc/ (double) s )
                    threshold[counter++] = ts*s*(d+1) - s*delta + 2*Pc -1;
                
                load[counter-1] = threshold[counter-1] / ((double) (t*d));
                
            }
            else
            {
                if ((Pc % MIN(t,d))==0)
                    delta = Pc / MIN(t,d);
                else
                    delta = ((int) (Pc / MIN(t,d)))+1;
                
                ss = s + delta;
                ds = d + delta;
                
                if (t==d &&  ((double) delta > Pc/ (double) s))
                    threshold[counter++] = ss*(t*t+1) -3;
                else
                    threshold[counter++] = t*d*ss + ss -1;
                
                load[counter-1] = threshold[counter-1] / ((double) (t*d));
            }
        }

    delete [] threshold;
    delete [] load;
    
    char filnavn[100];
    sprintf(filnavn,"Load_versus_threshold_P%d_Pc%d_m%d_n%d.txt",P,Pc,m,n);
    FILE *fp = fopen(filnavn,"w");
    
    fprintf(fp,"load threshold\n");
    for (int i=0;i<counter;i++)
        fprintf(fp,"%f %f\n",load[i],threshold[i]);
    fprintf(fp,"\n");
        
    fclose(fp);
}
