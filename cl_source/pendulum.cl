#define __OPENCL_C_VERSION__ 300

#define G __global
#define GE __generic
#define K __kernel void
typedef unsigned char byte;
typedef double real;

typedef struct _data{
	real p0,p1;//angle
	real v0,v1;//velo
	real m,l;//m0/m1,l1/l0
	real gx,gy;//gravity
	real e;//energy
}data;
const real t=0.02;
const real Pi=3.14159265358979323846246;

real hypernorm(real x){
	return x+sqrt(x*x+1);
}
inline real kinetic(GE data *p){
	real ct=cos(p->p0-p->p1);
	real c=p->v0,d=p->v1,e=p->m,f=p->l;
	return e*c*c+c*c+d*d*f*f+2*c*d*f*ct;
}
inline real potential(GE data *p){
	real s0=sin(p->p0),c0=cos(p->p0),s1=sin(p->p1),c1=cos(p->p1);
	real e=p->m,f=p->l;
	real gx=p->gx,gy=p->gy;
	return -2*(e*c0+gx*(s0+s1*f)+gy*(c0+c1*f));
}
// ca.m=7.0/9;
// ca.l=2.0/3;
K init(
	int sx,int sy,int w0,int w1,
	G real *px,G real *py,G data *a,G data *s)
{
	const int y=get_global_id(0),p=y*sx;
	for(int i=0;i<sx;++i){
		data ca=*a;
		real GE *g[8]={&ca.p0,&ca.p1,&ca.v0,&ca.v1,&ca.m,&ca.l,&ca.gx,&ca.gy};
		*g[w0]=px[i];
		*g[w1]=py[y];
		ca.m=hypernorm(-ca.m);
		ca.l=hypernorm(ca.l);
		ca.e=kinetic(&ca)+potential(&ca);
		s[p+i]=ca;
	}
}
typedef struct _slope{
	real a0,a1,v0,v1;
}slope;
slope emu(GE data *p,GE slope *s,real t){
	real a=p->p0,b=p->p1,c=p->v0,d=p->v1,e=p->m,f=p->l;
	if(t!=0){
		a+=s->v0*t;
		b+=s->v1*t;
		c+=s->a0*t;
		d+=s->a1*t;
	}
	real st=sin(a-b),ct=cos(a-b),s0=sin(a),s1=sin(b),c1=cos(b);
	real gx=p->gx,gy=p->gy;
	real a0=s0+st*(c*c*ct+f*d*d-s0*st+gx*s1+gy*c1)/(e+st*st);
	real a1=(c*c*st+a0*ct+gx*c1-gy*s1)/f;
	return (slope){-a0,a1,c,d};
}
inline void emulate(GE data *p){
	slope rk[4];
	rk[0]=emu(p,0,0);
	rk[1]=emu(p,rk,t/2);
	rk[2]=emu(p,rk+1,t/2);
	rk[3]=emu(p,rk+2,t);
	p->p0+=t*(rk[0].v0+(rk[1].v0+rk[2].v0)*2+rk[3].v0)/6;
	p->p1+=t*(rk[0].v1+(rk[1].v1+rk[2].v1)*2+rk[3].v1)/6;
	p->v0+=t*(rk[0].a0+(rk[1].a0+rk[2].a0)*2+rk[3].a0)/6;
	p->v1+=t*(rk[0].a1+(rk[1].a1+rk[2].a1)*2+rk[3].a1)/6;
	
	real k=kinetic(p);
	if(k>0.001){
		real rate=sqrt((p->e-potential(p))/k);
		p->v0*=rate;
		p->v1*=rate;
	}
	if(p->p0>=-100&&p->p0<=100){
		while(p->p0<-Pi)p->p0+=Pi+Pi;
		while(p->p0> Pi)p->p0-=Pi+Pi;
	}
	if(p->p1>=-100&&p->p1<=100){
		while(p->p1<-Pi)p->p1+=Pi+Pi;
		while(p->p1> Pi)p->p1-=Pi+Pi;
	}
}
inline void emulate2(GE data *p){
	const static real g=10;
	real a=p->p0,b=p->p1,c=p->v0,d=p->v1;
	real st=sin(a-b),ct=cos(a-b),s0=sin(a),s1=sin(b);
	real de=16-9*ct*ct;
	real v0=(2*c-3*ct*d)*6/de;
	real v1=(8*d-3*ct*c)*6/de;
	real a0=-(3*g*s0+v0*v1*st)/2;
	real a1=-(g*s1-v0*v1*st)/2;
	p->p0+=v0*t;
	p->p1+=v1*t;
	p->v0+=a0*t;
	p->v1+=a1*t;
}
K calc(int sz,int step,G data *s){
	const int id=get_global_id(0),p=id*sz;
	for(int i=p;i<p+sz;++i)
	for(int j=0;j<step;++j)emulate(s+i);
}
K norm(int sz,G data *s,G byte *tar){
	const int id=get_global_id(0),p=id*sz;
	for(int i=p;i<p+sz;++i){
		real x=s[i].p0,y=s[i].p1;
		int g=i*3;
		tar[g  ]=cos(x)*120+128;
		tar[g+1]=(sin(x)+cos(y))*60+128;
		tar[g+2]=sin(y)*120+128;
	}
}