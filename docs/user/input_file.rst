User manual
===========

.. image:: ../logo/blue.png
   :alt: REIMS illustration
   :width: 360

Input file description
----------------------

  REIMS input file is used to describe the model which will be calculated
  by REIMS software.
  
  Small example of input file is shown below:

  .. code:: yaml
    
    simulation:         
      simulation_end: 1000        # Simulate 1000s
      implicit_tolerance: 0.0005  # Tolerance for implicit solver
    
    write_results:            
      file: reims_output.h5       # Write results to file: reims_output.h5
    
    # Main model description is a list of components and their connections
    components:         # One "state" component and one "link" component
      - type: channel   # Type of the component - In this case channel
        id: pipe        # Name of the component has to be unique
        nodes: 200      # Number of computation cells
        length: 140.0   # Total length in m
        diameter: 10e-3 # Channel diameter 12mm
        initial: {p: 5.0e5, t: 4.3} # initial conditions P = 5bar and T = 4.3K
      - type: pump   # Type of the component - In this case pump
        m0: 2.0e-3   # Mass flow rate: 2 g/s
        link:        # 2 links: link 1 - pump inlet, link 2 - pump outlet
         - id: pipe  # inlet of the pump connected to outlet of the 'pipe'
           node: out # outlet pipe
         - id: pipe  # outlet of the pump connected to inlet of the 'pipe'
           node: in  # inlet pipe
  
  It describes 2 components:
  
    1. ``channel`` which is state component, and name (id) ``pipe``
    2. ``pump`` which is link component which links 2 ends of ``pipe``

  connected in a closed loop.

.. _top-level-configuration:

YAML config file structure
--------------------------

.. jsonschema:: ../../tools/reims_schema.json
   :lift_description: True
   :lift_definitions: True
   :auto_reference: True
   :auto_target: True
   :hide_key: /**/default

.. _friction-correlations-details:

Friction correlations details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Three use cases are supported:

**Case 1 — Built-in correlation with default parameters**

No ``friction_correlations`` section is needed. Reference the correlation
directly in the channel component:

.. code:: yaml

  components:
    - type: channel
      friction: blasius

Built-in names: ``blasius``, ``katheder``, ``central_spiral``.

``blasius``: f = alpha * Re\ :sup:`-beta` — Defaults: alpha=0.079, beta=0.25.

``katheder``: f = 0.25 * (19.5 / Re\ :sup:`beta` + alpha) / VoidFr\ :sup:`0.742` — Defaults: alpha=0.0231, beta=0.7953, VoidFr=0.297, Re_min=1000.
Returns zero for Re ≤ ``Re_min``; only laminar friction applies below this threshold. Set ``Re_min: 0`` to activate Katheder for all Re > 0.

``central_spiral``: f = 0.25 * alpha / Re\ :sup:`beta` — Defaults: alpha=0.36, beta=0.038.

where Re is the Reynolds number.

**Case 2 — Built-in correlation with custom parameters**

Define a named entry in ``friction_correlations`` using ``base`` to select
the built-in, then override the parameters:

.. code:: yaml

  friction_correlations:
    - friction: my_blasius
      base: blasius
      alpha: 0.04
      beta: 0.22

  components:
    - type: channel
      friction: my_blasius

**Case 3 — User-defined correlation via external DLL**

Provide a DLL that exports ``init_friction_ext`` and the correlation function.
Declare the library in ``external_libs`` and define the correlation name and
its parameters in ``friction_correlations``:

.. code:: yaml

  external_libs:
  - my_lib_friction

  friction_correlations:
    - friction: friction_test1
      alpha_test: 0.75

  components:
    - type: channel
      friction: friction_test1

The parameters defined under ``friction_correlations`` are passed to the DLL
at startup. See the `DLL developer guide <../../../dll/_build/html/index.html>`_ for the required interface.


.. _nusselt-correlations-details:

Nusselt correlations details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Three use cases are supported:

**Case 1 — Built-in correlation with default parameters**

No ``nusselt_correlations`` section is needed. Reference the correlation
directly in the component:

.. code:: yaml

  components:
    - type: fluidlink
      nusselt: pipe

Built-in names: ``pipe``, ``DBG``.

``pipe``: Nu = a * Re\ :sup:`b` * Pra\ :sup:`c` — Defaults: a=0.023, b=0.8, c=0.4.

``DBG``: NuT = a * Re\ :sup:`b` * Pra\ :sup:`c` * (T2/T1)\ :sup:`d` → Nu = max(NuL, NuT) — Defaults: a=0.0259, b=0.8, c=0.4, d=-0.716, NuL=8.235.

Re is the Reynolds number, Pra the Prandtl number, T1/T2 the temperatures (K) of the
two components. For ``channel`` (``thermal: temp``), T2 is the imposed wall temperature.

**Case 2 — Built-in correlation with custom parameters**

Define a named entry in ``nusselt_correlations`` using ``base``, then override the parameters:

.. code:: yaml

  nusselt_correlations:
    - nusselt: my_pipe
      base: pipe
      a: 0.02
      b: 0.85

  components:
    - type: fluidlink
      nusselt: my_pipe

**Case 3 — User-defined correlation via external DLL**

Declare the library in ``external_libs``, add a named entry with any custom parameters:

.. code:: yaml

  external_libs:
  - my_lib_nusselt

  nusselt_correlations:
    - nusselt: nusselt_test1
      param_test: 1.5

  components:
    - type: fluidlink
      nusselt: nusselt_test1

The parameters are passed to ``init_nusselt_ext`` at startup.
See the `DLL developer guide <../../../dll/_build/html/index.html>`_ for the required interface.


.. _materials-details:

Materials details
~~~~~~~~~~~~~~~~~~

Optional list of named materials available to ``strand``, ``solid``,
``solidlink``, and ``mesh2D`` components.
Three use cases are supported:

**Case 1 — Built-in material with default parameters**

No ``materials`` section is needed. Reference the built-in name directly
in the component:

.. code:: yaml

  components:
    - type: solid
      material: stainless_steel

Available built-in materials and their implemented properties:

- ``copper``:

  - density: 8960 kg/m³
  - resistivity

    .. container:: ref

      J. Simon, E.S. Drexler, R.P. Reed,
      Properties of Copper and Copper Alloys at Cryogenic Temperatures,
      NIST Monograph 177, Washington DC, 1992 (draft 1987 version including B dependence)

  - thermal conductivity

    .. container:: ref

      J. Simon, E.S. Drexler, R.P. Reed,
      Properties of Copper and Copper Alloys at Cryogenic Temperatures,
      NIST Monograph 177, Washington DC, 1992.
      B dependence added via a magnetoresistive term alpha\*B/T/L0
      (alpha~5e-11 Ohm.m/T, L0=2.44e-8 V2/K2) complementing the NIST dataset.

  - heat capacity

    .. container:: ref

      L. Dresner, Stability of Superconductors, Plenum Press, NY, 1995

- ``nb3sn``:

  - density: 8040 kg/m³
  - thermal conductivity

    .. container:: ref

      Fit of MATPRO data:
      L. Rossi, M. Sorbi, MATPRO: A Computer Library of Material Property at Cryogenic Temperature,
      INFN/TC-06/02, CARE-Note-2005-018-HHH.
      Original data from: H. Brechna, Superconducting Magnet Systems, Springer, 1973, p. 434.

  - heat capacity

    .. container:: ref

      ITER DRG1 Annex, Superconducting Material Database, Article 5,
      N 11 FDR 42 01-07-05 R 0.1, Thermal, Electrical and Mechanical
      Properties of Materials at Cryogenic Temperatures (internal ITER report).
      Original data from:
      V.D. Arp, Stability and Thermal Quenches in Force-Cooled Superconducting Cables,
      Superconducting MHD Magnet Design Conf., MIT, pp 142-157, 1980;
      G.S. Knapp, S.D. Bader, Z. Fisk, Phonon properties of A-15 superconductors
      obtained from heat capacity measurements, Phys. Rev. B, 13(9), pp 3783-3789, 1976.

  - strain

    .. container:: ref

      Linear model with a constant term and an electromagnetic contribution proportional to I x B.

  - critical temperature

    .. container:: ref

      L. Bottura, B. Bordini, Jc(B,T,epsilon) Parameterization for the ITER Nb3Sn Production,
      IEEE Trans. Appl. Sup., 19(2), 1477-1480, 2009

  - critical field

    .. container:: ref

      L. Bottura, B. Bordini, Jc(B,T,epsilon) Parameterization for the ITER Nb3Sn Production,
      IEEE Trans. Appl. Sup., 19(2), 1477-1480, 2009

  - critical current density

    .. container:: ref

      L. Bottura, B. Bordini, Jc(B,T,epsilon) Parameterization for the ITER Nb3Sn Production,
      IEEE Trans. Appl. Sup., 19(2), 1477-1480, 2009

  - current sharing temperature

- ``nbti``:

  - density: 6000 kg/m³
  - thermal conductivity

    .. container:: ref

      6th-order polynomial fit from MATPRO dataset:
      L. Rossi, M. Sorbi, MATPRO: A Computer Library of Material Property at Cryogenic Temperature,
      INFN/TC-06/02, CARE-Note-2005-018-HHH.
      Original data from: H. Brechna, Superconducting Magnet Systems, Springer, 1973, p. 424.

  - heat capacity

    .. container:: ref

      Elrod S.A., Miller J.R., Dresner L.,
      The specific heat of NbTi from 0-7T between 4.2 and 20K,
      Advances in Cryogenic Engineering Materials, Vol. 28, 1981.
      Extended to the full temperature range by smooth transition
      to a high-temperature asymptote of 400 J/kg/K.

  - strain

    .. container:: ref

      Linear model with a constant term and an electromagnetic contribution proportional to I x B.

  - critical temperature

    .. container:: ref

      M.S. Lubell, Scaling formulas for critical current and critical field for commercial NbTi,
      IEEE Trans. Mag., 19, 1983

  - critical field

    .. container:: ref

      M.S. Lubell, Scaling formulas for critical current and critical field for commercial NbTi,
      IEEE Trans. Mag., 19, 1983

  - critical current density

    .. container:: ref

      M.A. Green, Calculating the Jc, B, T Surface for Niobium Titanium Using a Reduced State Model,
      IEEE Trans. Mag., 25(2), 1989;
      G. Morgan, A Comparison of Two Analytic Forms for the Jc(B,T) surface, SSC-MD-218, 1989;
      L. Bottura, B. Bordini, Jc(B,T,epsilon) Parameterization for the ITER Nb3Sn Production,
      IEEE Trans. Appl. Sup., 19(2), 1477-1480, 2009;
      L. Zani, J.P. Serries, H. Cloez, M. Tena, E. Mossang,
      Jc(B,T) characterization of NbTi strands used in the ITER PF (Poloidal Field Coil)-relevant
      Insert and Full-scale sample, INIS-FR--2832, 2004

  - current sharing temperature

- ``stainless_steel``:

  - density: 7900 kg/m³
  - thermal conductivity

    .. container:: ref

      J.M. Poncet, CEA-Grenoble, EFDA CRYOLA task.

  - heat capacity

    .. container:: ref

      ITER DRG1 Annex, Superconducting Material Database, Article 5, N 11 FDR 42 01-07-05 R 0.1,
      Thermal, Electrical and Mechanical Properties of Materials at Cryogenic Temperatures (internal ITER report).
      Debye fit; for the model formula see for instance:
      L. Dresner, Stability of Superconductors, Plenum Press, NY, 1995.
      Original data from:
      J.M. Corsan and N.I. Mitchem, The Specific Heat of fifteen stainless steels in the temperature range 4K-30K,
      Cryogenics 19, p11-p16;
      J.M. Corsan and N.I. Mitchem, The Specific Heat of stainless steels between 4K and 300K,
      Proc. of the 6th ICEC, Grenoble, 1976;
      Aerospace Structural Metals Handbook, Metals and Ceramics Information Center,
      Battelle's Columbus Laboratories, Columbus, OH.

- ``glass_epoxy``:

  - density: 1948 kg/m³
  - thermal conductivity

    .. container:: ref

      Cubic polynomial fit of data from:
      M.B. Kasen, G.R. MacDonald, D.H. Beekman Jr and R.E. Schrmm,
      Mechanical, Electrical and Thermal Characterisation of G-10CR and G-11CR
      Glass-Cloth/Epoxy Laminates Between Room Temperature and 4K,
      National Bureau of Standards, Boulder, Colorado.

  - heat capacity

    .. container:: ref

      Power-law fit of data from:
      G. Hartwig, Low Temperature Properties of Resins and Their Correlations,
      Adv. Cryog. Eng., Vol. 24, 1978.

- ``glass_kapton_glass``:

  - density: 1800 kg/m³
  - thermal conductivity

    .. container:: ref

      J.M. Poncet, J.P. Arnaud, P. Saint Bonnet,
      Thermal conductivity of materials or sandwiches used for magnet insulation of ITER project,
      SBT report CT 12-42, February 2013.

  - heat capacity

    .. container:: ref

      Power-law fit of data from:
      G. Hartwig, Low Temperature Properties of Resins and Their Correlations,
      Adv. Cryog. Eng., Vol. 24, 1978.

**Built-in material parameters**

.. list-table::
   :widths: 15 85
   :header-rows: 0

   * - ``density``
     - Material density (kg/m³). Overrides the built-in default for this entry.
   * - ``RRR``
     - Residual Resistivity Ratio. Magnetoresistance correction factor in the copper resistivity model. Default: 100. *copper only*
   * - ``E0``
     - Electric field criterion (V/m). Reference field defining critical current in the power law. Default: 1.0e-5 V/m. *superconductors only*
   * - ``nPow``
     - Power law exponent. Exponent n in E = E0*(J/Jc)^n. Default: 5. *superconductors only*
   * - ``Bc20m``
     - Upper critical field at 0 K and zero intrinsic strain (T). Default: 29.39 T. *nb3sn only*
   * - ``Tc0m``
     - Critical temperature at zero field and zero intrinsic strain (K). Default: 16.48 K. *nb3sn only*
   * - ``nu``
     - Shape exponent for the temperature dependence of Bc2. nb3sn: 1.52, nbti: 1.7.
   * - ``Ca1``
     - First strain fitting constant. Default: 45.74. *nb3sn only*
   * - ``Ca2``
     - Second strain fitting constant. Default: 4.431. *nb3sn only*
   * - ``e0a``
     - Residual strain component. Default: 0.00232. *nb3sn only*
   * - ``emax``
     - Applied strain at which critical properties reach their maximum. Default: 0.0. *nb3sn only*
   * - ``C0``
     - Overall Jc scaling constant (A.T/m²). Default: 8.0771e10. *nb3sn only*
   * - ``p``
     - Low-field flux pinning exponent. nb3sn: 0.556, nbti: 0.98.
   * - ``q``
     - High-field flux pinning exponent. nb3sn: 1.698, nbti: 0.98.
   * - ``str1``
     - Strain constant term. nb3sn: 0.0060942, nbti: 0.00742.
   * - ``str2``
     - Strain electromagnetic coefficient (1/(A.T)). nb3sn: 1.0777e-9, nbti: 1.301e-9.
   * - ``Bc20``
     - Upper critical field at 0 K (T). Default: 13.72 T. *nbti only*
   * - ``Tc0_p``
     - Critical temperature at zero field (K). Default: 8.79 K. *nbti only*
   * - ``CC0``
     - Overall Jc scaling constant (A.T/m²). Default: 8.92534e11. *nbti only*
   * - ``n``
     - Exponent of the temperature-dependent factor (1-t^nu) in the Jc parameterization. Default: 1.96. *nbti only*

**Case 2 — Built-in material with custom parameters**

Define a named entry in ``materials`` using ``base`` to select the built-in,
then override only the parameters that differ from the defaults:

.. code:: yaml

  materials:
    - material: nb3sn_custom
      base: nb3sn
      Bc20m: 30.23
      Tc0m:  16.73
      E0:    1.0e-5
      nPow:  5

    - material: copper_rrr110
      base: copper
      RRR: 110.0

  components:
    - type: strand
      stabilizer:
        material: copper_rrr110
        area: 5.0e-7
      superconductor:
        material: nb3sn_custom
        area: 2.5e-7

**Case 3 — User-defined material via external DLL**

Provide a DLL that exports ``init_material_ext``.
Declare the library in ``external_libs`` and define the material name and
its parameters in ``materials``:

.. code:: yaml

  external_libs:
  - my_lib_material

  materials:
    - material: my_sc
      my_param: 1.0

  components:
    - type: strand
      superconductor:
        material: my_sc
        area: 2.5e-7

The parameters defined under ``materials`` are passed to the DLL at startup.
Only the properties actually implemented in the DLL are overridden — the
remaining properties fall back to the built-in values if ``base`` is also
provided, otherwise they remain unimplemented and trigger a runtime error if called.
See the `DLL developer guide <../../../dll/_build/html/index.html>`_ for the required interface.


YAML standard
-------------

To describe models YAML file is used.
Please see documentation of file format here:

https://yaml.org/

Supported standards are:

 * YAML 1.2.2
 * YAML 1.1
 * YAML 1.0

Extended with merge dictionary feature: ``<<:`` Which was not included in the
standard however it is very useful. Description in the draft for YAML 1.1:

https://yaml.org/type/merge.html

Fortunately many other system support it. So it became "semi" standard.


Include other YAML or JSON files
--------------------------------

Instead of putting any value in YAML file we replace it with dictionary contain:

.. code:: yaml

   include: file_name.json
  
The file will be included in this specific place. Files can be ``json`` or
``yaml``.

For example:

.. code:: yaml

   components:
   - type: channel
      id: pipe
      mesh: variable
      nodes: {include: nodes.yaml}

Please be aware that 3rd party software like language server or ``schema``
validator doesn't understand this feature and might issue an error however 
REIMS will work correctly.

Execution command lines
-----------------------

The following command lines allow the user to customize the configuration for parallel computing:

.. code:: yaml

    set OMP_DISPLAY_ENV=TRUE       # Shows OpenMP environment settings when program starts
    set OMP_PROC_BIND=close        # Binds threads near each other (to nearby cores)
    set OMP_PLACES=threads         # Places each OpenMP thread on a separate hardware thread
    set OMP_NUM_THREADS=96         # Uses 96 threads for OpenMP parallel regions - Number of threads to be adapted
    set MKL_NUM_THREADS=30         # Uses 30 threads for Intel MKL (Math Kernel Library) - Number of threads to be adapted
    set MKL_DEBUG_CPU_TYPE=5       # Simulates a specific CPU type (for debugging MKL)
    set MKL_ENABLE_INSTRUCTIONS=5  # Forces MKL to use a specific instruction set (like AVX-512)
    set MKL_DISPLAY_ENV=TRUE       # Shows MKL environment info at runtime

Here is the execution command line to apply in a command prompt, considering the executable reims.exe file as well as a given input.yaml file:

.. code:: yaml

    path_exe\reims.exe path_inp\input.yaml 
   
The keywords ``path_exe`` and ``path_inp`` refer to the paths of the executable file and the input file, respectively.

How to cite REIMS
-----------------

D. Furfaro, J. Kosek, A. Ovcharov, T. Schioler, R. Rotella, T. Luce,
A new fast and robust thermo-hydraulic code for ITER superconducting magnet simulation,
Cryogenics, Volume 144, 2024, 103978