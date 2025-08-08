# Architecture
- Top-level umbrella in workflow/Snakefile
- Stage modules in workflow/modules/<stage> (each has its own rules/envs/scripts)
- Long-lived references under resources/
- Outputs under results/<stage>/
