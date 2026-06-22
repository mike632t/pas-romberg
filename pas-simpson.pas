{
*  simpson.pas
* 
*  Copyright(C) 2026 - MT
*
*  Calculates  the definite integral for a function fn(x)  using  Simpson's
*  method. 
*
*  This  program is free software: you can redistribute it and/or modify it
*  under  the terms of the GNU General Public License as published  by  the
*  Free  Software Foundation, either version 3 of the License, or (at  your
*  option) any later version.
*
*  This  program  is distributed in the hope that it will  be  useful,  but
*  WITHOUT   ANY   WARRANTY;   without even   the   implied   warranty   of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
*  Public License for more details.
*
*  You  should have received a copy of the GNU General Public License along
*  with this program.  If not, see <http://www.gnu.org/licenses/>.
*
*  https://rosettacode.org/wiki/Numerical_integration/Romberg_integration
*
*  21 Jun 26   0.1   - Initial version - MT
* 
}
PROGRAM simpson(output);

CONST
   LOWER = -3.0;     { Lower bound of integration }
   UPPER = 3.0;      { Upper bound of integration }
   MAX   = 12;       { Maximum number of refinements }
   LIMIT = 5E-12;    { Convergence tolerance }

VAR
   a  : DOUBLE;      { Lower limit }
   b  : DOUBLE;      { Upper limit }
   r  : DOUBLE;      { Result }

FUNCTION Fn(x : DOUBLE) : DOUBLE;  { Function to be integrated }
BEGIN
   {Fn := sin(x);    { f(x) = sin(x) }
   {Fn := 1.0 / x;   { f(x) = 1 / x }
   Fn := exp(x);     {F(x) := exp(x)}
END;

FUNCTION simpson(FUNCTION op(x : DOUBLE) : DOUBLE; a, b : DOUBLE; max : INTEGER) : DOUBLE;

VAR
   h  : DOUBLE;      { Step size }
   s0 : DOUBLE;      { Odd-indexed sum }
   s1 : DOUBLE;      { Even-indexed sum }
   p  : DOUBLE;      { Previous estimate }
   c  : DOUBLE;      { Current estimate }
   d  : DOUBLE;      { Difference }
   i, j, n: INTEGER; { Loop counters }

BEGIN
   n := 2;  { Start with the minimum valid number of intervals (must be even) }
   h := (b - a) / n;
   
   p := op(a) + op(b) + 4.0 * op(a + h);  { Initial Simpson approximation using 2 intervals }
   p := p * h / 3.0;

   i := 1;  { Iteration counter }

   REPEAT
      i := i + 1;

      n := n * 2;  { Double number of intervals each iteration }
      h := (b - a) / n;  { Recompute step size }

      s0 := 0.0;
      s1 := 0.0;

      FOR j := 1 TO n - 1 DO  { Evaluate function ateach interval }
      BEGIN
         IF (j MOD 2 = 0) THEN
            s1 := s1 + op(a + j * h)
         ELSE
            s0 := s0 + op(a + j * h);
      END;

      c := op(a) + op(b) + 4.0 * s0 + 2.0 * s1;  { Simpson's rule f(a) + f(b) + 4*odd + 2*even }
      c := c * h / 3.0;
      d := ABS(c - p);  { Delta from previous iteration }
      p := c;  { Update previous estimate }

      WRITELN('I=', i:2, '     R= ', c:22:15);
   UNTIL (i >= max) OR (d < LIMIT);

   WRITELN;

   IF (i < MAX) THEN
      WRITELN('Converged early at I=', i:2)
   ELSE
      WRITELN('Max iterations reached');

   simpson := c;
END;

BEGIN
   a := LOWER;
   b := UPPER;
   r := simpson(Fn, a ,b , MAX);
   WRITELN('Integral =  ', r:22:15);
END.
