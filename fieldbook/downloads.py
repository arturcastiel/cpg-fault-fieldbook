"""Opt-in, checksum-verified HTTPS downloads. Never called by normal lessons."""

import hashlib
import tempfile
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import urlopen


def fetch(url: str, destination: Path, sha256: str, *, max_bytes: int) -> Path:
    """Fetch a licensed file only when caller supplies URL, size cap and hash."""
    if urlparse(url).scheme != "https":
        raise ValueError("only HTTPS downloads are accepted")
    if len(sha256) != 64 or any(c not in "0123456789abcdef" for c in sha256.lower()):
        raise ValueError("expected a 64-character SHA-256 digest")
    if max_bytes <= 0:
        raise ValueError("max_bytes must be positive")
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256()
    total = 0
    temporary: Path | None = None
    try:
        with urlopen(url, timeout=30) as response, tempfile.NamedTemporaryFile(
            dir=destination.parent, prefix=".download-", delete=False
        ) as out:
            temporary = Path(out.name)
            if urlparse(response.geturl()).scheme != "https":
                raise ValueError("download redirected outside HTTPS")
            while block := response.read(1024 * 1024):
                total += len(block)
                if total > max_bytes:
                    raise ValueError(f"download exceeds {max_bytes} bytes")
                digest.update(block)
                out.write(block)
        if digest.hexdigest() != sha256.lower():
            raise ValueError("download SHA-256 mismatch; destination unchanged")
        temporary.replace(destination)
        return destination
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
