project('CellConstructor', 'fortran', 'c', 'python',
  version: '1.5.0',
  license: 'MIT',
  # Opciones comunes de compilación para Fortran.
  # -cpp es esencial para el preprocesador Fortran.
  # -O2 es una optimización general que tenías en setup.py.
  default_options: ['fortran_args=-cpp', 'fortran_args=-O2']
)

# --- Dependencias del sistema y de Python ---
# Importa el módulo Python de Meson para manejar las extensiones y la instalación de Python.
python = import('python')
python_installation = python.find_installation()

# Define las dependencias de Python que numpy.distutils solía manejar.
# Meson las buscará en el entorno Python actual.
numpy_dep = python_installation.dependency('numpy', required: true)
ase_dep = python_installation.dependency('ase', required: true)
scipy_dep = python_installation.dependency('scipy', required: true)

# Busca las librerías BLAS y LAPACK en el sistema.
# Meson puede encontrar implementaciones como OpenBLAS, BLIS, Intel MKL, o las genéricas.
lapack_dep = dependency('lapack', required: true)
blas_dep = dependency('blas', required: true)

# Busca la librería OpenMP si es necesaria para la paralelización.
openmp_dep = dependency('openmp', required: true)


# --- Definición de las Extensiones Fortran ---

# Extensión 'symph'
symph_sources = files(
  'FModules/symdynph_gq_new.f90', 'FModules/symm_base.f90',
  'FModules/sgam_ph.f90', 'FModules/invmat.f90', 'FModules/set_asr.f90',
  'FModules/error_handler.f90', 'FModules/io_global.f90',
  'FModules/flush_unit.f90', 'FModules/symvector.f90',
  'FModules/fc_supercell_from_dyn.f90', 'FModules/set_tau.f90',
  'FModules/cryst_to_car.f90', 'FModules/recips.f90',
  'FModules/q2qstar_out.f90', 'FModules/rotate_and_add_dyn.f90',
  'FModules/trntnsc.f90', 'FModules/star_q.f90', 'FModules/eqvect.f90',
  'FModules/symm_matrix.f90', 'FModules/from_matdyn.f90',
  'FModules/interp.f90', 'FModules/q_gen.f90', 'FModules/smallgq.f90',
  'FModules/symmetry_high_rank.f90', 'FModules/unwrap_tensors.f90',
  'FModules/get_latvec.f90', 'FModules/contract_two_phonon_propagator.f90',
  'FModules/get_q_grid_fast.f90', 'FModules/kind.f90',
  'FModules/constants.f90', 'FModules/eff_charge_interp.f90',
  'FModules/get_translations.f90', 'FModules/get_equivalent_atoms.f90'
)
# Crea un módulo de extensión Python usando los fuentes Fortran y las dependencias.
py_symph_ext = python_installation.extension_module(
  'symph', # Nombre del módulo resultante (ej. symph.cpython-310-x86_64-linux-gnu.so)
  symph_sources,
  # Argumentos de compilación específicos para estos fuentes Fortran.
  # -cpp ya se definió en default_options, pero se puede añadir aquí si es necesario override.
  # fortran_args: ['-cpp'],
  dependencies: [numpy_dep, lapack_dep, blas_dep]
)

# Extensión 'secondorder'
secondorder_sources = files(
  'FModules/second_order_centering.f90',
  'FModules/second_order_ASR.f90'
)
py_secondorder_ext = python_installation.extension_module(
  'secondorder',
  secondorder_sources,
  # fortran_args: ['-cpp'],
  dependencies: [numpy_dep, lapack_dep, blas_dep]
)

# Extensión 'thirdorder'
thirdorder_sources = files(
  'FModules/third_order_centering.f90',
  'FModules/third_order_ASR.f90',
  'FModules/third_order_interpol.f90',
  'FModules/third_order_dynbubble.f90'
)
py_thirdorder_ext = python_installation.extension_module(
  'thirdorder',
  thirdorder_sources,
  # fortran_args: ['-cpp'],
  dependencies: [numpy_dep, lapack_dep, blas_dep]
)

# Extensión 'thermal_conductivity'
cond_sources = files(
  'FModules/get_scattering_q_grid.f90',
  'FModules/third_order_cond.f90',
  'FModules/third_order_cond_centering.f90',
  'FModules/get_lf.f90'
)
py_cond_ext = python_installation.extension_module(
  'thermal_conductivity',
  cond_sources,
  # Aquí se añade explícitamente '-fopenmp' para esta extensión específica.
  # "-O2" ya está en default_options.
  fortran_args: ['-fopenmp'],
  dependencies: [numpy_dep, lapack_dep, blas_dep, openmp_dep]
)

# --- Definición de la Extensión C 'cc_linalg' ---
# Lógica condicional para el archivo wrapper basada en la versión de Python.
wrapper_file = ''
if python_installation.version().version_compare('<3.0')
  wrapper_file = 'CModules/wrapper.c'
else
  wrapper_file = 'CModules/wrapper3.c'
endif

cc_modules_sources = files('CModules/LinAlg.c', wrapper_file)
py_cc_linalg_ext = python_installation.extension_module(
  'cc_linalg',
  cc_modules_sources,
  # Asumo que cc_linalg también necesita numpy, si no, se puede quitar.
  dependencies: [numpy_dep]
)

# --- Instalación de los módulos de extensión ---
# Registra cada módulo de extensión para que Meson sepa que debe instalarlos
# en el directorio de sitio-packages de Python.
python_installation.install_extmod(py_symph_ext)
python_installation.install_extmod(py_secondorder_ext)
python_installation.install_extmod(py_thirdorder_ext)
python_installation.install_extmod(py_cond_ext)
python_installation.install_extmod(py_cc_linalg_ext)


# --- Instalación del paquete Python 'cellconstructor' ---
# Esto instala los archivos .py del directorio 'cellconstructor'
# en el subdirectorio 'cellconstructor' dentro del sitio-packages de Python.
python_installation.install_sources(
  python_installation.find_sources('cellconstructor/*.py'),
  subdir: 'cellconstructor'
)

# Instalación de datos del paquete (archivos .dat)
# Se instalan dentro del directorio 'cellconstructor/SymData' en el sitio-packages.
install_data('cellconstructor/SymData/gamma_symm_op.dat',
             'cellconstructor/SymData/space_groups_crystal.dat',
             install_dir: python_installation.get_install_dir() / 'cellconstructor' / 'SymData')

# --- Instalación de los scripts ejecutables ---
# Instala los scripts directamente en el directorio de scripts del entorno Python.
install_scripts(
  'scripts/symmetrize_dynmat.py',
  'scripts/cellconstructor_test.py',
  'scripts/view_scf_atoms.py'
)
