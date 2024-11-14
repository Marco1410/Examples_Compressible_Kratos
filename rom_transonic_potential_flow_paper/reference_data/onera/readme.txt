This data is from:

Experimental:
    V. Schmitt and F. Charpin. Pressure Distributions on the
    ONERA-M6-Wing at Transonic Mach Numbers, Experimental Data
    Base for Computer Program Assessment. Tech. rep. Report of the
    Fluid Dynamics Panel Working Group 04, AGARD AR 138, 1979
    https://www.sto.nato.int/publications/AGARD/AGARD-AR-138/AGARD-AR-138.pdf

FV Potential solution:
    @article{LYU2017951,
    title = {A fast and automatic full-potential finite volume solver on Cartesian grids for unconventional configurations},
    journal = {Chinese Journal of Aeronautics},
    volume = {30},
    number = {3},
    pages = {951-963},
    year = {2017},
    issn = {1000-9361},
    doi = {https://doi.org/10.1016/j.cja.2017.03.001},
    url = {https://www.sciencedirect.com/science/article/pii/S1000936117300730},
    author = {Fanxi LYU and Tianhang XIAO and Xiongqing YU},
    keywords = {Cartesian grids, CUT-cell, Finite volume method, Full potential equation, Grid adaptation, Kutta condition, Non-penetration condition},
    abstract = {To meet the requirements of fast and automatic computation of subsonic and transonic aerodynamics in aircraft conceptual design, a novel finite volume solver for full potential flows on adaptive Cartesian grids is developed in this paper. Cartesian grids with geometric adaptation are firstly generated automatically with boundary cells processed by cell-cutting and cell-merging algorithms. The nonlinear full potential equation is discretized by a finite volume scheme on these Cartesian grids and iteratively solved in an implicit fashion with a generalized minimum residual (GMRES) algorithm. During computation, solution-based mesh adaptation is also applied so as to capture flow features more accurately. An improved ghost-cell method is proposed to implement the non-penetration wall boundary condition where the velocity-potential of a ghost cell is modified by an analytic method instead. According to the characteristics of the Cartesian grids, the Kutta condition is applied by specially computing the gradients on Kutta-faces without directly assigning the potential jump to cells adjacent wake faces, which can significantly improve the solution converging speed. The feasibility and accuracy of the proposed method are validated by several typical cases of sub/transonic flows around an ONERA M6 wing, a DLR-F4 wing-body, and an unconventional figuration of a blended wing body (BWB). The validation cases demonstrate a fast convergence with fully automatic grid treatment and computation, and the results suggest its capacity in application for aircraft conceptual design.}
    }

FE Potential solution:
    @phdthesis{dissertation,
        author = {López Canalejo, Iñigo Pablo},
        title = {A finite-element transonic potential flow solver with an embedded wake approach for aircraft conceptual design},
        year = {2022},
        school = {Technische Universität München},
        pages = {192},
        language = {en},
        abstract = {This dissertation presents a novel embedded wake approach for the solution of the full-potential equation using finite elements. The ultimate goal of this embedded approach is to effectively perform aircraft aeroelastic optimization at the early stages of aircraft design, taking transonic effects into account. The proposed approach is validated and verified on unstructured meshes in two and three dimensions.},
        keywords = {},
        note = {},
        url = {https://mediatum.ub.tum.de/1633175},
    }