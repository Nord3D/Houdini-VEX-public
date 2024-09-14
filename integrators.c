function vector samplev(int opinput; vector pos)
{
    int pr;
    vector uv;
    xyzdist(opinput, pos, pr, uv);
    
    return primuv(opinput, "v", pr, uv);
}



/*********************************************************/
/*********************  Runge-Kutta  *********************/

/*****   From http://www.mymathlib.com/diffeq/...    *****/

float h = chf("h");

/* Runge-Kutta-1, Euler */
function vector rk1 (vector pos; float h)
{
    vector v;
    v = samplev(1, pos);
    
    return v;
}

/* Runge-Kutta-2, Heun */
function vector rk2 (vector pos; float h)
{
    vector v,v1,v2;
    v1 = samplev(1, pos);
    v2 = samplev(1, pos + h * v1);
    v  = (v1 + v2)/2;
    
    return v;
}

/* Runge-Kutta-2, Improved Euler */
function vector rk2ie (vector pos; float h)
{
    vector v,v1,v2;
    v1 = samplev(1, pos);
    v2 = samplev(1, pos + h * v1);
    v  = v2;
    
    return v;
}

/* Runge-Kutta-2, Ralston */
function vector rk2ral (vector pos; float h)
{
    vector v,v1,v2;
    float h23 = 2. * h / 3.;
    float h4 = h / 4.;
    v1 = samplev(1, pos);
    v2 = samplev(1, pos + v1 * h23);
    v = h4 * (v1 + 3. * v2);
       
    return v;
}


/* Runge-Kutta-3, V1 */
function vector rk31 (vector pos; float h)
{
    vector v,v1,v2,v3;
    v1 = samplev(1, pos);
    v2 = samplev(1, pos  +  .5 * h * v1);
    v3 = samplev(1, pos  -  h * v1  +  2. * h * v2);
    v  = (v1  +  4 * v2  +  v3) / 6;
    
    return v;
}

/* Runge-Kutta-3, V2 */
function vector rk32 (vector pos; float h)
{
    vector v,v1,v2,v3;
    v1 = samplev(1, pos);
    v2 = samplev(1, pos  +       h * v1 / 3.);
    v3 = samplev(1, pos  +  2. * h * v2 / 3.);
    v  = (v1  +  3. * v3) / 4;
    
    return v;
}


/* Runge-Kutta-4, Classical */
function vector rk4 (vector pos; float h)
{
    vector v,v1,v2,v3,v4;
    float h_2 = h * .5;
    v1 = samplev(1, pos);
    v2 = samplev(1, pos + h_2 * v1);
    v3 = samplev(1, pos + h_2 * v2);
    v4 = samplev(1, pos +   h * v3);
    v  = (v1 + 2*(v2 + v3) + v4)/6;
    
    return v;
}

/* Runge-Kutta-4, 3/8 */
function vector rk4_38 (vector pos; float h)
{
    vector v,v1,v2,v3,v4;
    float h_2 = h * .5;
    v1 = samplev(1, pos);
    v2 = samplev(1, pos + h_2 * v1);
    v3 = samplev(1, pos + h_2 * v2);
    v4 = samplev(1, pos +   h * v3);
    v  = (v1 + 3*(v2 + v3) + v4)/8;
    
    return v;
}



/*********************************************************/
/*********************   B F E C C   *********************/

function vector bfecc1 (vector pos; float h)
{
    vector v, v1,v2,v3, p1,p2,p3;

    v1 = samplev(1, pos);
    p1 = pos  +  h * v1;
    
    v2 = samplev(1, p1);
    p2 = p1  -  h * v2;
    
    p3 = pos  +  (pos  -  p2) * 0.5;
    
    v3 = samplev(1, p3);

    v  = (p3  +  h * v3  -  pos) / h ;
    
    return v;
}

function vector bfecc_rk2 (vector pos; float h)
{
    vector v, v1,v2,v3, p1,p2,p3;

    v1 = rk2(pos, h);
    p1 = pos  +  h * v1;
    
    v2 = rk2(p1, -h);   /// Внимание, нельзя просто инвертировать результат RK, нужно в RK -h
    p2 = p1  -  h * v2;
    
    p3 = pos  +  (pos  -  p2) * 0.5;
    
    v3 = rk2(p3, h);

    v  = ((p3  +  h * v3)  -  pos) / h;
    
    return v;
}

function vector bfecc_rk4 (vector pos; float h)
{
    vector v, v1,v2,v3, p1,p2,p3;

    v1 = rk4(pos, h);
    p1 = pos  +  h * v1;
    
    v2 = rk4(p1, -h);   /// Внимание, нельзя просто инвертировать результат RK, нужно в RK -h
    p2 = p1  -  h * v2;
    
    p3 = pos  +  (pos  -  p2) * 0.5;
    
    v3 = rk4(p3, h);

    v  = ((p3  +  h * v3)  -  pos) / h;
    
    return v;
}

function vector bfecc_rk4_38 (vector pos; float h)
{   /***** rk4_38 - Better RK4 variant *****/
    vector v, v1,v2,v3, p1,p2,p3;

    v1 = rk4_38(pos, h);
    p1 = pos  +  h * v1;
    
    v2 = rk4_38(p1, -h);   /// Внимание, нельзя просто инвертировать результат RK, нужно в RK -h
    p2 = p1  -  h * v2;
    
    p3 = pos  +  (pos  -  p2) * 0.667;  /// Гиперкомпенсация ))
    
    v3 = rk4_38(p3, h);

    v  = ((p3  +  h * v3)  -  pos) / h;
    
    return v;
}

/* Идея использовать Runge-Kutta в качестве "строительных болков" для BFECC, отсюда:
*  https://www.shadertoy.com/view/wt33z2 
*/



/* TEST */

v@v_ = h * rk1(@P, h);
v@v_ = h * rk2(@P, h);
// v@v_ = h * rk2ie(@P, h);
// v@v_ = h * rk2ral(@P, h);
v@v_ = h * rk31(@P, h);
// v@v_ = h * rk32(@P, h);
v@v_ = h * rk4_38(@P, h);

v@v_ = h * bfecc1(@P, h);
v@v_ = h * bfecc_rk2(@P, h);
v@v_ = h * bfecc_rk4(@P, h);
v@v_ = h * bfecc_rk4_38(@P, h);
