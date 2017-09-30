#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
using namespace std;
int main(){
    ofstream file_out;
    float g = 1;
    float dt=.00001;
    int n;
   double t,t_step;
   cin>>n>>t>>t_step;
   cout<<t_step<<endl;
   double arr[n][13];
   for(int i=0;i<n;i++){
    double m,x,y,z,vx,vy,vz;
    cin>>m;
    cin>>x;
    cin>>y;
    cin>>z;
    cin>>vx;
    cin>>vy;
    cin>>vz;
    arr[i][0]=m;
    arr[i][1]=x;
    arr[i][2]=y;
    arr[i][3]=z;
    arr[i][4]=vx;
    arr[i][5]=vy;
    arr[i][6]=vz;
   }

if(t_step>0){
   for(double i=0;i<=t;i=i+t_step){
    for(int j=0;j<n;j++){
        for (int k=j+1;k<n;k++){
        double dvx1,dvx2,dvy1,dvy2,dvz1,dvz2,dx1,dx2,dy1,dy2,dz1,dz2;
        double r=sqrt((arr[j][1]-arr[k][1])*(arr[j][1]-arr[k][1])+(arr[j][2]-arr[k][2])*(arr[j][2]-arr[k][2])+(arr[j][3]-arr[k][3])*(arr[j][3]-arr[k][3]));
        double force=g*arr[j][0]*arr[k][0]/r/r;
        double dir[3]={(arr[j][1]-arr[k][1])/r*force,(arr[j][2]-arr[k][2])/r*force,(arr[j][3]-arr[k][3])/r*force};
        dvx1=-dir[0]*dt/arr[j][0];
        dvy1=-dir[1]*dt/arr[j][0];
        dvz1=-dir[2]*dt/arr[j][0];

        dvx2=dir[0]*dt/arr[k][0];
        dvy2=dir[1]*dt/arr[k][0];
        dvz2=dir[2]*dt/arr[k][0];

        dx1=arr[j][4]*dt;
        dy1=arr[j][5]*dt;
        dz1=arr[j][6]*dt;

        dx2=arr[k][4]*dt;
        dy2=arr[k][5]*dt;
        dz2=arr[k][6]*dt;

        arr[j][1]=arr[j][1]+dx1;
        arr[j][2]=arr[j][2]+dy1;
        arr[j][3]=arr[j][3]+dz1;
        arr[j][4]=arr[j][4]+dvx1;
        arr[j][5]=arr[j][5]+dvy1;
        arr[j][6]=arr[j][6]+dvz1;

        arr[k][1]=arr[k][1]+dx2;
        arr[k][2]=arr[k][2]+dy2;
        arr[k][3]=arr[k][3]+dz2;
        arr[k][4]=arr[k][4]+dvx2;
        arr[k][5]=arr[k][5]+dvy2;
        arr[k][6]=arr[k][6]+dvz2;//Calculating forces as pairs
            }
        }
        for(int f1=0;f1<n;f1++){
        cout<<arr[f1][1]<<" "<<arr[f1][2]<<" "<<arr[f1][3]<<endl;
        file_out.open("output.txt",ios_base::app);
        file_out<<arr[f1][1]<<" "<<arr[f1][2]<<" "<<arr[f1][3]<<" "<<endl;
        file_out.close();
        }}

   return 0;



}}
