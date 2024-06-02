# utl-sympy-exact-pdf-and-cdf-for-the-correlation-coefficient-given-bivariate-normals
Use the excact probabilty density function for the correlation coefficient to calaculate    the p-value for a random sample of size 30 and a realized value of 0.36983 for the          coorrelation coef.                                                                         
    %let pgm=utl-sympy-exact-pdf-and-cdf-for-the-correlation-coefficient-given-bivariate-normals;

            Use the excact probabilty density function for the correlation coefficient to calaculate
            the p-value for a random sample of size 30 and a realized value of 0.36983 for the
            coorrelation coef.

    github
    https://tinyurl.com/23s3hu87
    https://github.com/rogerjdeangelis/utl-sympy-exact-pdf-and-cdf-for-the-correlation-coefficient-given-bivariate-normals



    Related
    Repos on end

    /***********************************************************************************************************************************/
    /*                                        |                                                   |                                    */
    /*                                        |                                                   |                                    */
    /*             INPUT                      |                 PROCESS                           |          OUTPUT                    */
    /*                                        |                                                   |                                    */
    /* correlation coef = 0.36983 and n=30    | %utl_pybegin;                                     | P-val=0.0443(SAS FisherZ 0.0437)   */
    /*                                        | parmcards4;                                       |                                    */
    /* pdf =                                  | import sympy as sp                                | cdf = 2*r*gamma(n/2 - 1/2)         */
    /*   (2*sp.gamma((n-1)/2)                 | # Define the symbols                              | *hyper((1/2, 2 - n/2), (3/2,)      */
    /*   *sp.sqrt(1-r**2)**(n-4))             | n=sp.symbols('n',positive=True,integer=True)      | ,r**2*exp_polar(2*I*pi))           */
    /*   /(sp.sqrt(sp.pi)*sp.gamma((n-2)/2))  | r=sp.symbols('r',real=True) #Correlation          | /(sqrt(pi)*gamma(n/2 - 1))         */
    /*                                        | # DefinethePDFsymbolically                        | + 2*gamma(n/2 - 1/2)               */
    /* proc corr data=sashelp.cars(obs=30)    | pdf=(2*sp.gamma((n-1)/2)*sp.sqrt(1-r**2)**(n-4)) \| *hyper((1/2, 2 - n/2), (3/2,), 1)  */
    /*       fisher(biasadj=no);;             | /(sp.sqrt(sp.pi)*sp.gamma((n-2)/2))               | /(sqrt(pi)*gamma(n/2 - 1))         */
    /* var msrp;                              | #Integrate thePDF ove range[-1,r]forcdf           |                                    */
    /* with weight;                           | cdf=sp.integrate(pdf,(r,-1,r))                    |                                    */
    /* run;quit;                              | print(cdf)                                        |                                    */
    /*                                        | integral=sp.integrate(pdf,(r,.36983,1))           |                                    */
    /* I don't believe the fisher z transform | result=integral.subs({n:30}).evalf()              |                                    */
    /* is exact but it handles issue with     | print(result)                                     |                                    */
    /* small sample sizes                     | ;;;;                                              |                                    */
    /* that tend to be non-normal             | %utl_pyend;                                       |                                    */
    /*                                        |                                                   |                                    */
    /*      With                     p Value  |                                                   |                                    */
    /* Var  Var     N   Correlation  H0:Rho=0 |                                                   |                                    */
    /*                                        |                                                   |                                    */
    /* SRP  WEIGHT  30  0.36983       0.0437  |                                                   |                                    */
    /*                                        |                                                   |                                    */
    /***********************************************************************************************************************************/

    /*                   _
    (_)_ __  _ __  _   _| |_
    | | `_ \| `_ \| | | | __|
    | | | | | |_) | |_| | |_
    |_|_| |_| .__/ \__,_|\__|
            |_|
    */

    /***********************************************************************************************************************************/
    /*                                                                                                                                 */
    /*                                                                                                                                 */
    /*             INPUT                                                                                                               */
    /*                                                                                                                                 */
    /* correlation coef = 0.36983 and n=30                                                                                             */
    /*                                                                                                                                 */
    /* pdf =                                                                                                                           */
    /*   (2*sp.gamma((n-1)/2)                                                                                                          */
    /*   *sp.sqrt(1-r**2)**(n-4))                                                                                                      */
    /*   /(sp.sqrt(sp.pi)*sp.gamma((n-2)/2))                                                                                           */
    /*                                                                                                                                 */
    /* proc corr data=sashelp.cars(obs=30)                                                                                             */
    /*       fisher(biasadj=no);;                                                                                                      */
    /* var msrp;                                                                                                                       */
    /* with weight;                                                                                                                    */
    /* run;quit;                                                                                                                       */
    /*                                                                                                                                 */
    /* I don't believe the fisher z transform                                                                                          */
    /* is exact but it handles issue with                                                                                              */
    /* small sample sizes                                                                                                              */
    /* that tend to be non-normal                                                                                                      */
    /*                                                                                                                                 */
    /*      With                     p Value                                                                                           */
    /* Var  Var     N   Correlation  H0:Rho=0                                                                                          */
    /*                                                                                                                                 */
    /* SRP  WEIGHT  30  0.36983       0.0437                                                                                           */
    /*                                                                                                                                 */
    /***********************************************************************************************************************************/

    pdf =
       (2*sp.gamma((n-1)/2)
       *sp.sqrt(1-r**2)**(n-4))
       /(sp.sqrt(sp.pi)*sp.gamma((n-2)/2))

    proc corr data=sashelp.cars(obs=30)
          fisher(biasadj=no);;
    var msrp;
    with weight;
    run;quit;

    /*
     _ __  _ __ ___   ___ ___  ___ ___
    | `_ \| `__/ _ \ / __/ _ \/ __/ __|
    | |_) | | | (_) | (_|  __/\__ \__ \
    | .__/|_|  \___/ \___\___||___/___/
    |_|
    */

    %utl_pybegin;
    parmcards4;
    import sympy as sp
    # Define the symbols
    n=sp.symbols('n',positive=True,integer=True)
    r=sp.symbols('r',real=True) #Correlation
    # DefinethePDFsymbolically
    pdf=(2*sp.gamma((n-1)/2)*sp.sqrt(1-r**2)**(n-4)) \
    /(sp.sqrt(sp.pi)*sp.gamma((n-2)/2))
    #Integrate thePDF ove range[-1,r]forcdf
    cdf=sp.integrate(pdf,(r,-1,r))
    print(cdf)
    integral=sp.integrate(pdf,(r,.36983,1))
    result=integral.subs({n:30}).evalf()
    print(result)
    ;;;;
    %utl_pyend;

    /*           _               _
      ___  _   _| |_ _ __  _   _| |_
     / _ \| | | | __| `_ \| | | | __|
    | (_) | |_| | |_| |_) | |_| | |_
     \___/ \__,_|\__| .__/ \__,_|\__|
                    |_|
    */

    /**************************************************************************************************************************/
    /*                                                                                                                        */
    /* P-value=0.0442698125897017                                                                                             */
    /*                                                                                                                        */
    /* cdf=gamma(n/2 - 1/2)*hyper((1/2, 2 - n/2), (3/2,), r**2*exp_polar(2*I*pi))/(sqrt(pi)*gamma(n/2 - 1))                   */
    /* + 2*gamma(n/2 - 1/2)*hyper((1/2, 2 - n/2), (3/2,), 1)/(sqrt(pi)*gamma(n/2 - 1))                                        */
    /*                                                                                                                        */
    /**************************************************************************************************************************/

    /*
     _ __ ___ _ __   ___  ___
    | `__/ _ \ `_ \ / _ \/ __|
    | | |  __/ |_) | (_) \__ \
    |_|  \___| .__/ \___/|___/
             |_|
    */
    REPO
    -----------------------------------------------------------------------------------------------------------------------------------
    https://github.com/rogerjdeangelis/utl-calculating-the-cube-root-of-minus-one-with-drop-down-to-python-symbolic-math-sympy
    https://github.com/rogerjdeangelis/utl-distance-between-a-point-and-curve-in-sql-and-wps-pythony-r-sympy
    https://github.com/rogerjdeangelis/utl-maximum-liklihood-regresssion-wps-python-sympy
    https://github.com/rogerjdeangelis/utl-python-sympy-projection-of-the-intersection-of-two-parabolic-surfaces-onto-the-xy-plane-AI
    https://github.com/rogerjdeangelis/utl-r-python-compute-the-area-between-two-curves-AI-sympy-trapezoid
    https://github.com/rogerjdeangelis/utl-solve-a-system-of-simutaneous-equations-r-python-sympy
    https://github.com/rogerjdeangelis/utl-symbolic-algebraic-simplification-of-a-polynomial-expressions-sympy
    https://github.com/rogerjdeangelis/utl-sympy-technique-for-symbolic-integration-of-bivariate-density-function

    /*              _
      ___ _ __   __| |
     / _ \ `_ \ / _` |
    |  __/ | | | (_| |
     \___|_| |_|\__,_|

    */
