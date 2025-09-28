#include"headers/clform.h"
#include<cmath>

typedef double real;
source p0;
void InitialFractal(){p0.Load("cl_source\\pendulum.cl");}
struct data{
	real p0,p1;//angle
	real v0,v1;//velo
	real m,l;//m0/m1,l1/l0
	real gx,gy;//gravity
	real e;//energy
};

const real Pi=acos(-1);

const real t=0.00000045;
real hypernorm(real x){
	return x+sqrt(x*x+1);
}
inline real kinetic(data *p){
	real ct=cos(p->p0-p->p1);
	real c=p->v0,d=p->v1,e=p->m,f=p->l;
	return e*c*c+c*c+d*d*f*f+2*c*d*f*ct;
}
inline real potential(data *p){
	real s0=sin(p->p0),c0=cos(p->p0),s1=sin(p->p1),c1=cos(p->p1);
	real e=p->m,f=p->l;
	real gx=p->gx,gy=p->gy;
	return -2*(e*c0+gx*(s0+s1*f)+gy*(c0+c1*f));
}
void initdata(data *p){
	p->m=hypernorm(-p->m);
	p->l=hypernorm(p->l);
	p->e=kinetic(p)+potential(p);
}
struct slope{
	real a0,a1,v0,v1;
};
slope emu(data *p,slope *s,real t){
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
inline void emulate(data *p){
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
inline void emulate2(data *p){
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

const real add[8]={Pi,Pi,Pi,Pi,3,3,3,3};
struct fractal{
	kernel init,calc,norm;
	int sx,sy,wx,wy;//w:enum direct
	int depth;
	data s,s0;//current,simulator
	real *p[8];
	real *px,*py;
	memory mc,ms,mx,my,mca;
	int sz,st;//pixel per item,step at cpu
	bool state;//0:multi at gpu 1:single at cpu
	void sstate(){
		state^=1;
	}
	void initial(int x,int y,int w0,int w1,int de,int size,int step,int step2){
		sx=x;sy=y;
		wx=w0;wy=w1;
		depth=de;
		sz=size;
		st=step2;
		state=0;
		p[0]=&s.p0;p[1]=&s.p1;p[2]=&s.v0;p[3]=&s.v1;
		p[4]=&s.m ;p[5]=&s.l ;p[6]=&s.gx;p[7]=&s.gy;
		for(int i=0;i<7;++i)*p[i]=0;*p[7]=1;
		init.Load(p0,"init");
		calc.Load(p0,"calc");
		norm.Load(p0,"norm");
		px=new real[sx];
		py=new real[sy];
		mc.Buffout(sizeof(data));
		ms.Buffout(sizeof(data)*sx*sy);
		mx.Buffout(sizeof(real)*sx);
		my.Buffout(sizeof(real)*sy);
		mca.Buffout(sizeof(byte)*sx*sy*3);
		init.Setint(0,sx);
		init.Setint(1,sy);
		init.Setint(2,wx);
		init.Setint(3,wy);
		init.Setptr(4,mx);
		init.Setptr(5,my);
		init.Setptr(6,mc);
		init.Setptr(7,ms);
		calc.Setint(0,sz);
		calc.Setint(1,step);
		calc.Setptr(2,ms);
		norm.Setint(0,sz);
		norm.Setptr(1,ms);
		norm.Setptr(2,mca);
	}
	void start(){
		real rate=1;
		for(int i=0;i<depth;++i)rate/=2;
		for(int i=0;i<sx;++i)px[i]=*p[wx]+(i-(sx>>1))*rate*add[wx];
		for(int i=0;i<sy;++i)py[i]=*p[wy]+(i-(sy>>1))*rate*add[wy];
		mc.Buffwrite(&s);
		mx.Buffwrite(px);
		my.Buffwrite(py);
		init.call1(0,sy,256);
		s0=s;
		initdata(&s0);
		puts("Reset");
	}
	void play(byte *buf){
		if(state){
			for(int i=0;i<st;++i)emulate(&s0);
		}
		else{
			calc.call1(0,sx*sy/sz,256,1);
			if(buf){
				norm.call1(0,sx*sy/sz,256);
				mca.Buffread(buf);
			}
		}
	}
	void offset(fractal *s,int x,int y){
		if(x<0||x>=s->sx||y<0||y>=s->sy)return;
		*p[s->wx]=s->px[x];
		*p[s->wy]=s->py[y];
		start();
	}
	void zoom(int a){
		if(depth+a<4)return;
		depth+=a;
		if(a)start();
	}
};