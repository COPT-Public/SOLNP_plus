/* Installation test: min (x-1)^2+(y-2)^2, x+y=3, -5<=x,y<=5.
 * Analytic optimum: (1,2), objective 0. No derivatives are supplied.
 */
#include "solnp_util.h"
#include <math.h>
#include <stdio.h>
static void cost(double *x, double *v, int n) {
    (void)n;
    v[0]=(x[0]-1)*(x[0]-1)+(x[1]-2)*(x[1]-2);
    v[1]=x[0]+x[1]-3;
}
int main(int argc, char **argv) {
    if(argc!=2) { fprintf(stderr,"Usage: test_nlp problem.txt\n"); return 2; }
    FILE *f=fopen(argv[1],"r");
    double x0,y0,lb,ub;
    if(!f || fscanf(f,"%lf %lf %lf %lf",&x0,&y0,&lb,&ub)!=4) return 2;
    fclose(f);
    srand(20260930);
    SOLNPSettings *s=calloc(1,sizeof(*s));
    SOLNPSol *sol=calloc(1,sizeof(*sol));
    SOLNPInfo *info=calloc(1,sizeof(*info));
    SOLNPProb *p=SOLNP(init_prob)(2,0,1);
    SOLNP(set_default_settings)(s,2);
    s->tol=1e-7; s->tol_con=1e-7; s->max_iter=100;
    s->maxfev=5000; s->qpsolver=1; s->rs=0;
    p->p0=malloc(2*sizeof(double)); p->p0[0]=x0; p->p0[1]=y0;
    p->pbl=malloc(2*sizeof(double)); p->pbu=malloc(2*sizeof(double));
    for(int i=0;i<2;i++) { p->pbl[i]=lb; p->pbu[i]=ub; }
    SOLNP_PLUS(p,s,sol,info,cost,NULL,NULL,NULL,NULL);
    double x=sol->p[0],y=sol->p[1];
    double obj=(x-1)*(x-1)+(y-2)*(y-2), res=fabs(x+y-3);
    int ok=(sol->status==1)&&isfinite(x)&&isfinite(y)&&isfinite(sol->obj)&&
        fabs(x-1)<1e-3&&fabs(y-2)<1e-3&&obj<1e-6&&res<1e-5&&
        fabs(sol->obj-obj)<1e-6&&x>=lb-1e-5&&x<=ub+1e-5&&y>=lb-1e-5&&y<=ub+1e-5;
    printf("COINOR_RESULT {\"status\":%d,\"x\":[%.17g,%.17g],\"objective\":%.17g,\"constraint_residual\":%.17g,\"pass\":%s}\n",sol->status,x,y,obj,res,ok?"true":"false");
    SOLNP(free_sol)(sol); SOLNP(free_info)(info); SOLNP(free_stgs)(s); SOLNP(free_prob)(p);
    return ok?0:1;
}
