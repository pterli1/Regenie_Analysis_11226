"""DNAnexus platform, for UK Biobank and other DNAnexus-hosted biobanks.

[`dx`][fedcat.platforms.dnanexus.dx] wraps dxpy;
[`DxPlatform`][fedcat.platforms.dnanexus.platform.DxPlatform] implements the platform
interface with it.
"""

from fedcat.platforms.dnanexus.platform import DxCfg, DxPlatform

__all__ = ["DxCfg", "DxPlatform"]
