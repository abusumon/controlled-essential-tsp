#include <bits/stdc++.h>
using namespace std;
int n; vector<double> X,Y;
double d(int a,int b){double dx=X[a]-X[b],dy=Y[a]-Y[b];return sqrt(dx*dx+dy*dy);}
double L(const vector<int>&t){double s=0;for(int k=0;k<n;k++)s+=d(t[k],t[(k+1)%n]);return s;}
bool tryMove(vector<int>&t,int m,int x,int y){
  if(x==y) return false;
  int i,j; if(m==2){i=x;j=y;} else {i=min(x,y);j=max(x,y);}
  int pI=(i+n-1)%n,nI=(i+1)%n,pJ=(j+n-1)%n,nJ=(j+1)%n;
  double g;
  if(m==0){ int ci=t[i],cj=t[j],pi=t[pI],ni=t[nI],pj=t[pJ],nj=t[nJ];
    if(j-i==1) g=d(pi,ci)+d(cj,nj)-d(pi,cj)-d(ci,nj);
    else if(i==0&&j==n-1) g=d(pj,cj)+d(ci,ni)-d(pj,ci)-d(cj,ni);
    else g=d(pi,ci)+d(ci,ni)+d(pj,cj)+d(cj,nj)-d(pi,cj)-d(cj,ni)-d(pj,ci)-d(ci,nj);
    if(g>1e-9){swap(t[i],t[j]);return true;} return false;}
  if(m==1){ if(j<i+2||(i==0&&j==n-1)) return false;
    int ci=t[i],ci1=t[nI],cj=t[j],cj1=t[nJ];
    g=d(ci,ci1)+d(cj,cj1)-d(ci,cj)-d(ci1,cj1);
    if(g>1e-9){reverse(t.begin()+i+1,t.begin()+j+1);return true;} return false;}
  int ci=t[i],pi=t[pI],ni=t[nI],cj=t[j];
  if(i<j){ if(i==0&&j==n-1) return false; int nj=t[nJ];
    g=d(pi,ci)+d(ci,ni)+d(cj,nj)-d(pi,ni)-d(cj,ci)-d(ci,nj);
    if(g>1e-9){for(int k=i;k<j;k++)t[k]=t[k+1];t[j]=ci;return true;} return false;}
  if(pJ==i) return false; int pj=t[pJ];
  g=d(pj,cj)+d(pi,ci)+d(ci,ni)-d(pj,ci)-d(ci,cj)-d(pi,ni);
  if(g>1e-9){for(int k=i;k>j;k--)t[k]=t[k-1];t[j]=ci;return true;} return false;
}
int main(int argc,char**argv){
  ifstream f(argv[1]); string line; while(getline(f,line)){istringstream s(line);double a,b,c;string extra; if(s>>a>>b>>c && !(s>>extra)){X.push_back(b);Y.push_back(c);}}
  n=X.size(); int starts=atoi(argv[2]);
  mt19937_64 rng(12345);
  vector<int> init(n);
  for(int st=0;st<starts;st++){
    iota(init.begin(),init.end(),0); shuffle(init.begin(),init.end(),rng);
    printf("%d,%d,random,%.3f\n",n,st,L(init));
    const char* nm[3]={"swap","2opt","ins"};
    for(int m=0;m<3;m++){
      vector<int> t=init; bool imp=true; vector<int> ord(n); iota(ord.begin(),ord.end(),0);
      while(imp){ imp=false; shuffle(ord.begin(),ord.end(),rng);
        for(int a:ord) for(int y=0;y<n;y++) if(tryMove(t,m,a,y)) imp=true; }
      printf("%d,%d,%s,%.3f\n",n,st,nm[m],L(t)); if(getenv("DUMP")){FILE*q=fopen((string(getenv("DUMP"))+"_"+nm[m]+"_"+to_string(st)+".txt").c_str(),"w");for(int c:t)fprintf(q,"%d\n",c);fclose(q);} fflush(stdout);
    }
  }
}
