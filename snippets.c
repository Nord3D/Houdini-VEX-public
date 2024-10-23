
/*** MATH CONSTANTS ***/
#define NM_PI           3.14159265358979323846
#define NM_1_3          0.33333333333333333333      /// 1/3
#define NM_2_3          0.66666666666666666667      /// 2/3
#define NM_SQRT3_4      0.43301270189221932338      /// sqrt(3) / 4
#define NM_3SQRT3_4     1.29903810567665797014      /// 3 * sqrt(3) / 4     <Area of a Unit Equilateral Triangle>
#define NM_SQRT2_12     0.11785113019775792073      /// sqrt(2) / 12        <Volume of a Regular Tetrahedron with Unit edge length>



vector2 rot90(vector2 v)
{
    vector2 r = v;
    r[0] = v[1];
    r[1] = -v[0];
    return r;
}

int linesintersection(vector2 p11, p12, p21, p22, p_inter)
    /// Compute two lines intersection point in 2D
    /// p11, p12 - first line; p21, p22 - second line points; p_inter - resultant intersection point
    /// From https://www.topcoder.com/community/competitive-programming/tutorials/geometry-concepts-line-intersection-and-its-applications/
{
    float A1 = p12.y - p11.y;
    float B1 = p11.x - p12.x;
    float C1 = A1 * p11.x + B1 * p11.y;
    
    float A2 = p22.y - p21.y;
    float B2 = p21.x - p22.x;
    float C2 = A2 * p21.x + B2 * p21.y;
    
    float det = A1 * B2 - A2 * B1;
    if (det != 0){
        p_inter.x = (B2 * C1  -  B1 * C2) / det;
        p_inter.y = (A1 * C2  -  A2 * C1) / det;
        return 1;   /// Intersection point exist
    }else{
        return -1;  /// Lines are parallel
    }
}

float skewlinesdist(vector P11, P12, P21, P22, C1, C2)
/// P1X, P2X - lines points
/// C1, C2 - nearest points on the 1- and 2- lines
/// From https://en.wikipedia.org/wiki/Skew_lines#Nearest_points
{
    vector D1 = normalize(P12 - P11);
    vector D2 = normalize(P22 - P21);
    
    vector N = cross(D1, D2);
    vector N1 = cross(D1, N);
    vector N2 = cross(D2, N);
    
    vector S = P21 - P11;
    float d1 = dot( S, N2) / dot(D1, N2);
    float d2 = dot(-S, N1) / dot(D2, N1);
    
    /*** R E T U R N ***/
    C1 = P11 + d1 * D1;
    C2 = P21 + d2 * D2;
    return distance(C1, C2);
}

int segmentsintersect(vector P11,P12, P21,P22, P)
/// P1X, P2X - segments points
/// P - intersection point
{
    vector D1 = normalize(P12 - P11);
    vector D2 = normalize(P22 - P21);
    
    vector N = cross(D1, D2);
    vector N1 = cross(D1, N);
    vector N2 = cross(D2, N);
    
    vector S = P21 - P11;
    float d1 = dot( S, N2) / dot(D1, N2);
        if ( d1 < 0 || d1 > length(P12 - P11) ) return 0;
    float d2 = dot(-S, N1) / dot(D2, N1);
        if ( d2 < 0 || d2 > length(P22 - P21) ) return 0;

    P = P11 + d1 * D1;
    return 1;
}

int rayplaneintersection(vector r, r0, pn, p0; float d)
    /// r - ray (unit vector), r0 - ray starting point, pn - plane normal, p0 - plane point 
    /// d - resultant ray distance from starting point to intersection point with plane
    /// From https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection#Algebraic_form
{
    float denominator = dot(r, pn);
    if (denominator == 0){
        return 0;                               /// parallel case
    }else{
        d = dot(p0 - r0, pn) / denominator;  
        return 1;                               /// Intersection point exist
    }
}


vector orthogonalize(vector N, T)
{
    vector B = cross(N, T);
    return normalize(cross(T, B));
}

vector rodrigues(vector 
                        v,  // Vector to rotate
                        v1, // Source vector, normalized
                        v2  // Target vector, normalized
                    )
{
    vector axis = cross(v1, v2);
    vector axisn = normalize(axis);
    float cosq = dot(v1, v2);
    return v * cosq + cross(axis, v) + axisn * dot(axisn, v) * (1. - cosq);
}

vector rodihedral(vector v1,    // Vector to rotate
                         v2     // Target vector, normalized
                   )
{
    vector v1n = normalize(v1);
    vector axis = cross(v1n, v2);
    float cosq = dot(v1n, v2);
    return v1 * cosq + cross(axis, v1);
}

vector rodihedral(vector v1,    // Vector to rotate
                         v2;    // Target vector, normalized
                  float mix)
{
    vector v1n = normalize(v1);
    vector axis = cross(v1n, v2);
    float angle = acos(dot(v1n, v2)) * mix;
    axis = normalize(axis) * sin(angle);
    float cosq = cos(angle);
    return v1 * cosq + cross(axis, v1);
}




vector dirproject(
                  vector a; /// vector to project
                  vector b; /// vector to project onto, normalized
                  vector d; /// projection direction, normalized
                  )
/*https://math.stackexchange.com/questions/3027898/projection-of-vector-into-an-axis-along-a-direction*/
{
    /* Remove from a it's projection onto d */
    a -= d * dot(a, d);
    /* Then return a, projected on b orthogonal to a */
    return b * dot(a, a) / dot(a, b);
}




/****************************************************************
* Все эти циклы делают одно и то-же:                            *
* создают 5 новых точек со сдвигом по Y                         *
****************************************************************/
vector P = 0;
int n = 5;

/* for_multistring */
for ( int i = 0; i < n; ++i ) {
    P.y += .1;
    addpoint(0, P);
}
/* for_singlestring */
for ( int i = 0; i < n; ++i ) addpoint(0, P += {0,.1,0});
/* В цикле for безразлично, префиксный или постфиксный 
 * инкремент/декремент используется.
 */

/* while_multistring */
while ( counter-- ) {
    P.y += .1;
    addpoint(0, P);
}
/* while_singlestring */
while ( counter-- ) addpoint(0, P += {0,.1,0});
/* В циклах while и do-while, важно, префиксный или постфиксный 
 * инкремент/декремент используется.
 * В данном примере использован постфиксный декремент, поскольку
 * он возвращает предыдущее значение. 
 */

/* do_while_multistring */
do {
    P.y += .1;
    addpoint(0, P);
} while ( --counter );
/* do_while_singlestring */
do addpoint(0, P += {0,.1,0}); while ( --counter );
/* В циклах while и do-while, важно, префиксный или постфиксный 
 * инкремент/декремент используется.
 * В данном примере использован префиксный декремент, поскольку
 * он возвращает обновлённое (уменьшенное на 1) значение. 
 */


