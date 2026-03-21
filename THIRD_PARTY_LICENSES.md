ROCCO is distributed under the MIT license in [LICENSE](LICENSE).

This distribution vendors source from [HTSlib](https://www.htslib.org/) and
includes a native counting backend adapted from the
[Consenrich](https://github.com/nolan-h-hamilton/Consenrich) package.

The following third-party components are included:

- HTSlib
  - source path: `vendor/htslib`
  - upstream license file: `vendor/htslib/LICENSE`
  - packaged license copy: `HTSLIB_LICENSE.txt`
  - license summary: files outside `cram/` are under the MIT/Expat license and
    files in `cram/` are under the modified 3-clause BSD license

- htscodecs
  - source path: `vendor/htslib/htscodecs`
  - upstream license file: `vendor/htslib/htscodecs/LICENSE.md`
  - packaged license copy: `HTSCODECS_LICENSE.md`
  - license summary: BSD-style license with some files noted as public domain
    or CC0-derived in the upstream license text

- Consenrich adaptations
  - source path: `rocco/native/ccounts_backend.c`, `rocco/native/ccounts_backend.h`, `rocco/native/baseline_backend.c`, `rocco/native/baseline_backend.h`, `rocco/_baseline.c`, `rocco/inference.py`
  - upstream project license file: [Consenrich LICENSE](https://github.com/nolan-h-hamilton/Consenrich/blob/main/LICENSE)
  - license summary: MIT
  - notes: ROCCO adapts the native counting backend and the cross-fit Whittaker local baseline implementation from Consenrich

Please note:

- ROCCO keeps its own MIT license
- Bundled third-party code keeps its original licenses
- Source and binary redistributions should preserve upstream copyright notices,
  license terms, and disclaimers
