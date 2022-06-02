# B-spline Rationnelles Non uniformes et rendu

# sommaire

- [introduction](#introduction)

- [principe du B-spline](#principe-du-b-spline)
 
- [technique d'implementation](#technique-dimplementation)
 
- [Comparaison](#comparaison)
 
- [performance](#performance)



# introduction

La spline de base rationnelle non uniforme ( NURBS ) est un modèle mathématique utilisant des splines de base (B-splines) qui est couramment utilisé en infographie pour représenter des courbes et des surfaces . Il offre une grande flexibilité et précision pour gérer à la fois les formes analytiques (définies par des formules mathématiques courantes ) et modélisées 





# principe du B-spline





# technique d'implementation
- # Knot()
## Knot c'est la fonction qui genere le knot Vector contenant des valeurs qui seront prises par t selon des intervales donnes

c            = order of the basis function

n            = the number of defining polygon vertices

nplus2       = index of x() for the first occurence of the maximum knot vector value

nplusc       = maximum value of the knot vector -- $n + c$

x()          = array containing the knot vector


```c
void knot(int n, int c, int x[])
{
  int nplusc,nplus2,i;

  nplusc = n + c;
  nplus2 = n + 2;

  x[1] = 0;
  for (i = 2; i <= nplusc; i++){
  if ( (i > c) && (i < nplus2) )
      x[i] = x[i-1] + 1;
  else
     x[i] = x[i-1];
  }
}

```

- # basis()

## l'implementation de l'algorithme de "cox boor" pour generer les valeurs des fonctions de base
 c        = order of the B-spline basis function
 
 d        = first term of the basis function recursion relation
 
 e        = second term of the basis function recursion relation
 
 npts     = number of defining polygon vertices
 
 n[]      = array containing the basis functions
 
 n[1] contains the basis function associated with B1 etc.
       
 nplusc   = constant -- npts + c -- maximum number of knot values
 
 t        = parameter value
 
 temp[]   = temporary array
 
 x[]      = knot vector
 
 ```c
 void basis(int c, float t, int npts, int x[], float n[])
{
    int nplusc;
    int i, k;
    float d, e;
    float temp[20000];
    nplusc = npts + c;
    /* calculate the first order basis functions n[i][1]    */
    for (i = 1; i <= nplusc - 1; i++)
    {
        if ((t >= x[i]) && (t < x[i + 1]))
            temp[i] = 1;
        else
            temp[i] = 0;
    }
    /* calculate the higher order basis functions */
    for (k = 2; k <= c; k++)
    {
        for (i = 1; i <= nplusc - k; i++)
        {
            if (temp[i] != 0) /* if the lower order basis function is zero skip the                         calculation */
                d = ((t - x[i]) * temp[i]) / (x[i + k - 1] - x[i]);
            else
                d = 0;
            if (temp[i + 1] != 0) /* if the lower order basis function is zero skip the calculation  */
                e = ((x[i + k] - t) * temp[i + 1]) / (x[i + k] - x[i + 1]);
            else
                e = 0;

            temp[i] = d + e;
        }
    }

    if (t == (float)x[nplusc])
    { /*    pick up last point    */
        temp[npts] = 1;
    }

    /* put in n array   */

    for (i = 1; i <= npts; i++)
    {
        n[i] = temp[i];
    }
}
 ```
 
 

- # bspline()

## Fonction utilisant les deux fonctions precedentes pour generer des valeurs des coordonnees des points de la courbe Bspline finale 
b[]        = array containing the defining polygon vertices

b[1] contains the x-component of the vertex
        
b[2] contains the y-component of the vertex
        
b[3] contains the z-component of the vertex
        
k           = order of the \bsp basis function

nbasis      = array containing the basis functions for a single value of t

nplusc      = number of knot values

npts        = number of defining polygon vertices

p[,]        = array containing the curve points

 p[1] contains the x-component of the point
        
p[2] contains the y-component of the point
        

p[3] contains the z-component of the point
        
p1          = number of points to be calculated on the curve

t           = parameter value 0 <= t <= 1

x[]         = array containing the knot vector
```c++
void bspline(int npts, int k, int p1, float b[], float p[])
{
    int i, j, icount, jcount;
    int i1;
    int x[25000] = {0}; /* allows for 20 data points with basis function of order 5 */
    int nplusc;
    float step;
    float t;
    float nbasis[25000] = {0.0};
    float temp;
    nplusc = npts + k;
    /*  zero and re-dimension the knot vector and the basis array */
    /* generate the uniform open knot vector */
    knot(npts, k, x);
    icount = 0;
    /*calculate the points on the bspline curve */
    t = 0;
    step = ((float)x[nplusc]) / ((float)(p1 - 1));
    for (i1 = 1; i1 <= p1; i1++)
    {
        if ((float)x[nplusc] - t < 5e-6)
        {
            t = (float)x[nplusc];
        }
        basis(k, t, npts, x, nbasis); /* generate the basis function for this value of t */
        for (j = 1; j <= 3; j++)
        { /* generate a point on the curve */
            jcount = j;
            p[icount + j] = 0.;
            for (i = 1; i <= npts; i++)
            { /* Do local matrix multiplication */
                temp = nbasis[i] * b[jcount];
                for (int count = 0; count < 12; count++)
                    count++;
                p[icount + j] = p[icount + j] + temp;
                jcount = jcount + 3;
            }
        }
        icount = icount + 3;
        t = t + step;
    }
}
```

# Comparaison





# performance
