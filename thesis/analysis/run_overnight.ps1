<#
    run_overnight.ps1 -- unattended driver for the remaining thesis experiments.

    Runs, in order:
      a. remaining spectral null seeds (1,2,3)
      b. node2vec seed-variance baseline: Generic and Liver FC, two new seeds each
      c. node2vec null seeds (2,3)
      d. Phase 3 predict-all inference on one spectral null seed

    Behaviour:
      * Every step SKIPS work whose output already exists and is complete
        (>= 6,000 trials for LOO; non-empty ranks file for inference).
      * Every step writes its own timestamped log under thesis/analysis/out/logs/.
      * A failing step is logged and the driver CONTINUES to the next step.
      * Uses the project's .venv312 interpreter. Nothing in src/, src_prediction/,
        data/ or outputs/ is modified; all new output lands under thesis/analysis/out/.

    Examples (-Command form; quote the path and use the call operator &)
      powershell -ExecutionPolicy Bypass -Command "& 'C:\Users\chavy\OneDrive\python\Kinase_Substrate_HGN\thesis\analysis\run_overnight.ps1'"
      powershell -ExecutionPolicy Bypass -Command "& 'C:\Users\chavy\OneDrive\python\Kinase_Substrate_HGN\thesis\analysis\run_overnight.ps1' -Steps d"
      powershell -ExecutionPolicy Bypass -Command "& 'C:\Users\chavy\OneDrive\python\Kinase_Substrate_HGN\thesis\analysis\run_overnight.ps1' -Steps b,c"

    Measured runtimes (node2vec LOO is dominated by the trial loop, ~3.1-3.3 s
    per trial over 6,909 trials, i.e. ~6 h per run; spectral LOO ~62 min;
    Phase 3 null inference ~10-20 min).
#>

param(
    [string[]]$Steps = @('a', 'b', 'c', 'd')
)

$ErrorActionPreference = 'Continue'

# Resolve project root as the parent of thesis\analysis\
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$Root = Resolve-Path (Join-Path $ScriptDir '..\..')
Set-Location $Root

$Py = Join-Path $Root '.venv312\Scripts\python.exe'
if (-not (Test-Path $Py)) { Write-Host "FATAL: interpreter not found at $Py"; exit 1 }

$LogDir = Join-Path $Root 'thesis\analysis\out\logs'
New-Item -ItemType Directory -Force -Path $LogDir | Out-Null
$NullDir = Join-Path $Root 'thesis\analysis\out\null'
$SeedDir = Join-Path $Root 'thesis\analysis\out\seed_variance'
$InfDir  = Join-Path $Root 'thesis\analysis\out\null_inference'

$Master = Join-Path $LogDir ('overnight_{0}.log' -f (Get-Date -Format 'yyyyMMdd_HHmmss'))

function Say($msg) {
    $line = "[{0}] {1}" -f (Get-Date -Format 'yyyy-MM-dd HH:mm:ss'), $msg
    Write-Host $line
    Add-Content -Path $Master -Value $line -Encoding utf8
}

function Test-Complete($kind, $path, $minTrials) {
    if ($minTrials) { & $Py -m thesis.analysis.is_complete $kind $path $minTrials 2>$null }
    else            { & $Py -m thesis.analysis.is_complete $kind $path 2>$null }
    return ($LASTEXITCODE -eq 0)
}

function Invoke-Step($name, $checkKind, $checkPath, $moduleArgs) {
    if (Test-Complete $checkKind $checkPath 6000) {
        Say "SKIP  $name  (output already complete)"
        return
    }
    $log = Join-Path $LogDir ("{0}_{1}.log" -f $name, (Get-Date -Format 'yyyyMMdd_HHmmss'))
    Say "START $name  -> $log"
    $sw = [Diagnostics.Stopwatch]::StartNew()
    try {
        & $Py @moduleArgs *>&1 | Tee-Object -FilePath $log
        $code = $LASTEXITCODE
    } catch {
        $code = 1
        Add-Content -Path $log -Value ("EXCEPTION: " + $_.Exception.Message) -Encoding utf8
    }
    $sw.Stop()
    $mins = [math]::Round($sw.Elapsed.TotalMinutes, 1)
    if ($code -eq 0) { Say "DONE  $name  ($mins min)" }
    else             { Say "FAIL  $name  (exit $code, $mins min) -- continuing" }
}

Say "=== overnight driver starting. steps requested: $($Steps -join ',') ==="
Say "python: $Py"

# ---------------------------------------------------------------- step a
if ($Steps -contains 'a') {
    Say "--- STEP A: spectral null seeds ---"
    foreach ($s in 1, 2, 3) {
        Invoke-Step "a_null_spectral_seed$s" 'loo' (Join-Path $NullDir "null_spectral_seed${s}_results.csv.gz") `
            @('-m', 'thesis.analysis.null_network', '--embedding', 'spectral', '--seeds', "$s")
    }
}

# ---------------------------------------------------------------- step b
if ($Steps -contains 'b') {
    Say "--- STEP B: node2vec seed-variance baseline ---"
    # Order: one seed across BOTH networks before the second seed, so a partial
    # night still yields a complete generic-vs-liver_fc pair to compare.
    foreach ($pair in @(@('generic', 7), @('liver_fc', 7), @('generic', 13), @('liver_fc', 13))) {
        $net = $pair[0]; $s = $pair[1]
        Invoke-Step "b_seedvar_${net}_seed$s" 'loo' (Join-Path $SeedDir "seedvar_${net}_node2vec_seed${s}_results.csv.gz") `
            @('-m', 'thesis.analysis.seed_variance', '--network', $net, '--seed', "$s")
    }
}

# ---------------------------------------------------------------- step c
if ($Steps -contains 'c') {
    Say "--- STEP C: node2vec null seeds ---"
    # Seed 1 is included: the original seed-1 run used batch_words=4 and was
    # quarantined under out/null/_INVALID_batch_words4/. It must be redone with
    # batch_words=10000 to be comparable to the real Liver FC node2vec run.
    foreach ($s in 1, 2, 3) {
        Invoke-Step "c_null_node2vec_seed$s" 'loo' (Join-Path $NullDir "null_node2vec_seed${s}_results.csv.gz") `
            @('-m', 'thesis.analysis.null_network', '--embedding', 'node2vec', '--seeds', "$s")
    }
}

# ---------------------------------------------------------------- step d
if ($Steps -contains 'd') {
    Say "--- STEP D: Phase 3 inference on spectral null seed 1 ---"
    Invoke-Step 'd_null_inference_spectral_seed1' 'ranks' (Join-Path $InfDir 'nullinfer_spectral_seed1_ranks.csv.gz') `
        @('-m', 'thesis.analysis.null_inference', '--seed', '1', '--embedding', 'spectral')
}

# ---------------------------------------------------------------- wrap up
Say "--- regenerating null comparison tables ---"
try { & $Py -m thesis.analysis.null_compare *>&1 | Tee-Object -FilePath (Join-Path $LogDir 'null_compare.log') | Out-Null }
catch { Say "null_compare failed: $($_.Exception.Message)" }

Say "=== overnight driver finished ==="
Say "Master log: $Master"
