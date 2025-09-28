#include"headers/gdiform.h"
#include"pendulum.h"

int sx,sy,size,step,step2;
window w[4];
fractal f[4];
rgb cross,back,r0,r1;
bool pause[4];

void getdata(){
	sx=1024;sy=768;
	size=16;step=12;
	step2=200000;
	cross={0,0,0};
	back={0,0,0};
	r0={255,64,96};
	r1={32,128,255};
}
void keys(int id,u64 wp){
	if(wp=='Q'||wp=='q')pause[id]^=1;
	if(wp=='E'||wp=='e')f[id].sstate();
	if(wp==VK_SPACE)w[(id+1)%4].ready(1);
	int a=0;
	if(wp=='W'||wp=='w')a=1;
	if(wp=='S'||wp=='s')a=-1;
	f[id].zoom(a);
	if(wp=='R'||wp=='r'){
		*f[id].p[f[id].wx]=0;
		*f[id].p[f[id].wy]=0;
		f[id].start();
	}
}
void click(int id,int x,int y,bool right){
	if(right){
		for(int i=0;i<4;++i)
			if(i!=id)f[i].offset(&f[id],x,y);
	}
	else f[id].offset(&f[id],x,y);
}
void show(window *g,fractal *f){
	int dx=sx>>1,dy=sy>>1;
	if(f->state){
		for(int j=0;j<sy;++j)
		for(int i=0;i<sx;++i)g->drawpoint(i,j,back);
		real r=dx<dy?dx:dy;
		r=r*5/6;
		real l0=r/(1+f->s0.l),l1=r-l0;
		real x0=l0*sin(f->s0.p0),y0=l0*cos(f->s0.p0),
			 x1=x0+l1*sin(f->s0.p1),y1=y0+l1*cos(f->s0.p1);
		g->drawline({(real)dx,(real)dy},{x0+dx,y0+dy},r0);
		g->drawline({x0+dx,y0+dy},{x1+dx,y1+dy},r1);
	}
	else{
		for(int i=-4;i<=4;++i)
			g->drawpoint(i+dx,dy,cross);
		for(int i=-4;i<=4;++i)
			g->drawpoint(dx,i+dy,cross);
	}
	g->paste();
}
void play(){
	bool flag=1;
	for(int i=0;i<4;++i)if(!pause[i]){
		if(flag)puts("Next frame:");
		f[i].play(w[i].data);
		show(&w[i],&f[i]);
		flag=0;
	}
	if(flag)Sleep(30);
}
int main(){
	GetCLAPI();
	InitialCL();
	InitialFractal();
	getdata();
	f[0].initial(sx,sy,0,1,8,size,step,step2);
	f[1].initial(sx,sy,2,3,8,size,step,step2);
	f[2].initial(sx,sy,4,5,8,size,step,step2);
	f[3].initial(sx,sy,6,7,8,size,step,step2);
	for(int i=0;i<4;++i)f[i].start();
	w[0].InitialGDI(100,30 ,sx,sy,"Plane in 6D - 0");
	w[1].InitialGDI(200,90 ,sx,sy,"Plane in 6D - 1");
	w[2].InitialGDI(300,150,sx,sy,"Plane in 6D - 2");
	w[3].InitialGDI(400,210,sx,sy,"Plane in 6D - 3");
	for(int i=3;i>=0;--i)w[i].ready(1);
	pause[0]=0;
	pause[1]=pause[2]=pause[3]=1;
	play();
	while(w[0].end&&w[1].end&&w[2].end&&w[3].end){
		message();
		if(w[0].fc||w[1].fc||w[2].fc||w[3].fc)play();
		else Sleep(30);
	}
	return 0;
}