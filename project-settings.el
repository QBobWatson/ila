(setq-local projectile-project-type 'scons)
(setq-local projectile-project-compilation-cmd "nix develop -c scons")
(setq-local compilation-read-command nil)

(unless (boundp 'ila-mode)
  (defvar ila-mode-map (make-sparse-keymap))
  (define-key ila-mode-map (kbd "C-c C-c") 'projectile-compile-project)
  (define-minor-mode ila-mode
    "Minor mode for this project."
    :keymap ila-mode-map))

(when (eq major-mode 'nxml-mode)
  (ila-mode 1))
