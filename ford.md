---
src_dir: src
output_dir: docs/fpm-ford
project: odepack
summary: ODEPACK - A Systematized Collection of ODE Solvers
project_github: https://github.com/jacobwilliams/odepack
project_download:
author: Alan C. Hindmarsh as modified by John S. Urban, Jacob Williams
author_email: urbanjost@comcast.net
github: https://github.com/jacobwilliams/odepack
media_dir: docs/images
exclude_dir: archive
             FODDER
             example
css: docs/local.css
display: public
         protected
source: true
extensions: f90
            inc
proc_internals: true
sort: permission-alpha
graph: true
favicon: docs/images/favicon.ico
print_creation_date: true
extra_mods: iso_fortran_env:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html
            tomlf:https://toml-f.github.io/toml-f
creation_date: %Y-%m-%d %H:%M %z
md_extensions: markdown.extensions.toc
               markdown.extensions.smarty
---

{!README.md!}
