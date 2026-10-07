<#
  run_all.ps1 -- runs the analysis on Windows.

      .\run_all.ps1 -Setup         install the R packages into env\rlib, create a Python
                                   environment with Jupyter, register the R kernel
      .\run_all.ps1 -Execute       run notebooks 1, 2 and 3, in this order
      .\run_all.ps1 -Analysis      run notebook 2 only
      .\run_all.ps1 -Validation    run notebook 3 only
      .\run_all.ps1                setup, then notebooks 1, 2 and 3

  Requires R 4.4.3 (with Rtools44) and Python 3. The R code uses this file to find the
  project root, so it must stay in the root folder.
#>
param([switch]$Setup, [switch]$Execute, [switch]$Analysis, [switch]$Validation)
$ErrorActionPreference = "Stop"
$root = $PSScriptRoot
$env:RENAL_IRI_ROOT = $root
$Rscript = (Get-Command Rscript.exe -ErrorAction SilentlyContinue).Source
if (-not $Rscript) {
  $Rscript = @(
      "$env:LOCALAPPDATA\Programs\R\R-4.4.3\bin\Rscript.exe",
      "$env:ProgramFiles\R\R-4.4.3\bin\Rscript.exe",
      "${env:ProgramFiles(x86)}\R\R-4.4.3\bin\Rscript.exe"
    ) | Where-Object { Test-Path $_ } | Select-Object -First 1
  if (-not $Rscript) { throw "Rscript.exe not found. Install R 4.4.3 (+ Rtools44) or add it to PATH." }
}
$venv = "$root\env\py"
if ($Analysis -or $Validation) { $Execute = $true }
if (-not $Setup -and -not $Execute) { $Setup = $true; $Execute = $true }

if ($Setup) {
  Write-Host "== [1/3] Installing R packages into env\rlib ==" -ForegroundColor Cyan
  $env:R_LIBS_USER = "$root\env\rlib"
  & $Rscript "$root\R\00_install.R"

  Write-Host "== [2/3] Creating the Python environment with Jupyter ==" -ForegroundColor Cyan
  if (-not (Test-Path "$venv\Scripts\python.exe")) { py -3 -m venv $venv }
  & "$venv\Scripts\python.exe" -m pip install --upgrade pip setuptools wheel
  & "$venv\Scripts\python.exe" -m pip install jupyterlab nbconvert nbclient ipykernel jupyter-client nbformat

  Write-Host "== [3/3] Registering the R kernel ==" -ForegroundColor Cyan
  & $Rscript "$root\R\01_register_kernel.R"
}

if ($Execute) {
  $env:JUPYTER_DATA_DIR = "$root\env\jupyter"
  $env:JUPYTER_PATH     = "$venv\share\jupyter"
  $env:R_LIBS_USER      = "$root\env\rlib"
  Set-Location $root
  $books = if ($Validation) { @("03_validation.ipynb") }
           elseif ($Analysis) { @("02_analysis.ipynb") }
           else { @("01_datasets.ipynb", "02_analysis.ipynb", "03_validation.ipynb") }
  foreach ($nb in $books) {
    Write-Host "== Executing notebook\$nb ==" -ForegroundColor Cyan
    & "$venv\Scripts\jupyter.exe" nbconvert --to notebook --execute --inplace `
        --ExecutePreprocessor.timeout=36000 `
        --ExecutePreprocessor.kernel_name=ir-renaliri `
        "$root\notebook\$nb"
    if ($LASTEXITCODE -ne 0) { throw "notebook\$nb failed (exit $LASTEXITCODE)" }
  }
  Write-Host "Done. See results\*.csv and figures\*.png." -ForegroundColor Green
}
