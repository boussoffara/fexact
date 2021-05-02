//
// Created by Alaa on 01/05/2021.
//

int imax2(int x, int y)
{
    return (x < y) ? y : x;
}


int imin2(int x, int y)
{
    return (x < y) ? x : y;
}

double fmax2(double x, double y)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(y))
		return x + y;
#endif
    return (x < y) ? y : x;
}


double fmin2(double x, double y)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(y))
		return x + y;
#endif
    return (x < y) ? x : y;
}
