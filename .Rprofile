# activate renv
source("renv/activate.R")

# VS Code Session Watcher + httpgd (robust for tmux/Remote-SSH)
if (interactive()) {
  # tmux sometimes clobbers TERM_PROGRAM; fix it if we detect VS Code
  if (Sys.getenv("VSCODE_IPC_HOOK_CLI") != "" && Sys.getenv("TERM_PROGRAM") != "vscode") {
    Sys.setenv(TERM_PROGRAM = "vscode")
  }

  # Source the watcher with working dir set to the script’s folder
  init_path <- file.path(Sys.getenv("HOME"),
                         ".vscode-server/extensions/reditorsupport.r-2.8.6/R/session/init.R")
  if (file.exists(init_path)) {
    source(init_path, chdir = TRUE)  # this calls init_first(), which sets .First.sys = init_last
    # We are in startup, so R will call .First.sys() automatically a bit later.
  } else {
    message("VS Code R init script not found at: ", init_path)
  }

  # Start httpgd per-session so plots land in VS Code’s Plot viewer
  try(httpgd::hgd(silent = TRUE, port = 0), silent = TRUE)
}
