#include "functions.h"

using namespace std;

Cosine::Cosine(){                                  // Cosine empty constructor
   
   m_a = 1;
   m_b = 1;
   m_c = 0;

}

Cosine::Cosine(double a, double b, double c){      // Cosine constructor

   m_a = a;
   m_b = b;
   m_c = c;

}

Parabola::Parabola(){                              // Parabola empty constructor

   m_a = 1;
   m_b = 0;
   m_c = 0;

}

Parabola::Parabola(double a, double b, double c){  // Parabola constructor
   
   m_a = a;
   m_b = b;
   m_c = c;

}

SquareRoot::SquareRoot(){    // square function constructor
   
   m_a = 0;
   m_b = 1;
   m_c = 0;

}

SquareRoot::SquareRoot(double a, double b, double c){    // square function constructor
   
   m_a = a;
   m_b = b;
   m_c = c;

}