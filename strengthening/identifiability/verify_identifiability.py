"""Check archived-input integrity, saved-fit reproduction and inactive stages."""
import json
from pathlib import Path
import runpy
import shutil
import tempfile
from run_identifiability import discover_validation, sha256

root = Path(__file__).resolve().parent
validation = discover_validation()
result = json.loads((root / "identifiability_results.json").read_text())
checks = {}
checks["original_input_hashes_unchanged"] = all(
    sha256(validation / name) == digest
    for name, digest in result["provenance"]["input_sha256"].items())
checks["script_hash_matches_record"] = (
    sha256(root / "run_identifiability.py") == result["provenance"]["script_sha256"])
checks["calibration_subsets_only"] = (
    result["gitis"]["calibration_n"] == 9 and
    result["alem"]["calibration_n"] == 26 and
    result["held_out_responses_used"] is False)
checks["original_fits_not_replaced"] = result["original_fits_replaced"] is False
checks["gitis_saved_SSE_reproduced"] = abs(
    result["gitis"]["baseline_metrics"]["SSE"] -
    result["gitis"]["saved_objective_SSE_N160"]) < 1e-12
checks["alem_saved_MSE_reproduced"] = abs(
    result["alem"]["baseline_metrics"]["SSE"] / 26 -
    result["alem"]["saved_fit_mean_squared_error_N120"]) < 1e-12
inactive = result["gitis"]["inactive_threshold_screen"]
checks["gitis_inactive_threshold_predictions_unchanged"] = all(
    row["max_absolute_prediction_change"] < 1e-12 for row in inactive)
checks["gitis_threshold_variants_remain_inactive"] = all(
    row["maximum_B"] < row["B2"] < row["capacity"] and
    row["cap_events"] == 0 for row in inactive)
profiles = (result["gitis"]["B1_objective_profile"] +
            result["alem"]["capacity_objective_profile"] +
            result["alem"]["one_sided_beta_objective_profile"])
checks["all_local_profile_optimizers_succeeded"] = all(
    row["optimizer_success"] for row in profiles)
checks["no_profile_cap_event"] = all(row["cap_events"] == 0 for row in profiles)
checks["all_profile_states_admissible"] = all(
    row["minimum_carrier"] >= -1e-12 and
    row["maximum_balance_error"] < 1e-10 for row in profiles)
checks["plots_present"] = all(
    (root / ("identifiability_diagnostics." + suffix)).stat().st_size > 1000
    for suffix in ("png", "pdf"))
with tempfile.TemporaryDirectory(prefix="portable_check_", dir=root) as temporary:
    package = Path(temporary)
    copied_script = package / "strengthening/identifiability/run_identifiability.py"
    copied_script.parent.mkdir(parents=True)
    shutil.copyfile(root / "run_identifiability.py", copied_script)
    target_validation = package / "validation"
    (target_validation / "gitis").mkdir(parents=True)
    for name in ("validation_solver.py", "gitis/fit_gitis_results.json"):
        shutil.copyfile(validation / name, target_validation / name)
    namespace = runpy.run_path(str(copied_script), run_name="portable_path_check")
    checks["portable_archive_directory_discovery"] = (
        namespace["discover_validation"]() == target_validation)
record = {"checks": checks, "all_passed": all(checks.values())}
(root / "verification_checks.json").write_text(json.dumps(record, indent=2) + "\n")
print(json.dumps(record, indent=2))
if not record["all_passed"]:
    raise SystemExit("A diagnostic verification check failed; inspect the record")
