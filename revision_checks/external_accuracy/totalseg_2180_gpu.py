"""TotalSegmentator 2.18.0 (task total, default 1.5 mm models baked in) + DICOM/RTSTRUCT tooling."""

import modal

META = {
    "packages": ["TotalSegmentator", "torch", "pydicom", "SimpleITK", "nibabel", "scikit-image"],
    "gpu_default": "A10G",
    # job-time host: TCIA NBIA API (per-case CT + RTSTRUCT download)
    "egress_domains": ["services.cancerimagingarchive.net"],
}

_ENV = {
    "TOTALSEG_HOME_DIR": "/opt/totalseg",
    "PYTHONUNBUFFERED": "1",
}


def build(*, secrets=None):
    img = (
        modal.Image.debian_slim(python_version="3.11")
        .apt_install("libgl1", "libxrender1", "libglib2.0-0", "libxext6")
        .pip_install(
            "TotalSegmentator==2.18.0",
            "pydicom",
            "SimpleITK",
            "nibabel",
            "scikit-image",
            "scipy",
            "pandas",
            "requests",
        )
        .env(_ENV)
        .run_commands(
            # pre-download the default (non-fast) 'total' task weights at build time
            "totalseg_download_weights -t total",
            # disable outbound usage statistics (job sandbox is egress-fenced)
            "python -c \"from totalsegmentator.config import setup_totalseg, set_config_key; "
            "setup_totalseg(); set_config_key('send_usage_stats', False)\"",
            "python -c \"import importlib.metadata as m; print('TotalSegmentator', m.version('TotalSegmentator'), "
            "'torch', m.version('torch'), 'nnunetv2', m.version('nnunetv2'))\"",
            "pip freeze > /opt/pip_freeze.txt",
            "ls -R /opt/totalseg | head -50",
        )
    )
    return img, {}, dict(_ENV)
