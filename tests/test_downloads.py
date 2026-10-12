import hashlib
import io
import json

import pytest

from fieldbook.cli import main
from fieldbook.downloads import fetch


class Response(io.BytesIO):
    def __init__(self, data: bytes, url: str = "https://example.invalid/data"):
        super().__init__(data)
        self.url = url

    def geturl(self) -> str:
        return self.url


def test_checksum_verified_and_atomic(monkeypatch, tmp_path):
    data = b"example public data"
    monkeypatch.setattr("fieldbook.downloads.urlopen",
                        lambda url, timeout: Response(data))
    destination = tmp_path / "data.bin"
    digest = hashlib.sha256(data).hexdigest()
    assert fetch("https://example.invalid/data", destination, digest, max_bytes=100) == destination
    assert destination.read_bytes() == data
    with pytest.raises(ValueError, match="mismatch"):
        fetch("https://example.invalid/data", destination, "0" * 64, max_bytes=100)
    assert destination.read_bytes() == data


def test_refuses_unverified_or_oversized_download(monkeypatch, tmp_path):
    destination = tmp_path / "out"
    with pytest.raises(ValueError, match="HTTPS"):
        fetch("http://example.invalid/data", destination, "0" * 64, max_bytes=100)
    monkeypatch.setattr("fieldbook.downloads.urlopen",
                        lambda url, timeout: Response(b"too long"))
    with pytest.raises(ValueError, match="exceeds"):
        fetch("https://example.invalid/data", destination, "0" * 64, max_bytes=1)
    assert not destination.exists()


def test_redirect_must_stay_https(monkeypatch, tmp_path):
    monkeypatch.setattr(
        "fieldbook.downloads.urlopen",
        lambda url, timeout: Response(b"data", "http://example.invalid/plaintext"),
    )
    with pytest.raises(ValueError, match="redirected"):
        fetch("https://example.invalid/data", tmp_path / "out", "0" * 64, max_bytes=100)


def test_opt_in_cli_records_source_and_license(monkeypatch, tmp_path):
    data = b"redistributable fixture"
    monkeypatch.setattr("fieldbook.downloads.urlopen",
                        lambda url, timeout: Response(data))
    destination = tmp_path / "fixture"
    digest = hashlib.sha256(data).hexdigest()
    assert main(["fetch", "https://example.invalid/fixture", "--sha256", digest,
                 "--max-bytes", "100", "--license", "CC0",
                 "--attribution", "Example Publisher", "--out", str(destination)]) == 0
    provenance = json.loads((tmp_path / "fixture.source.json").read_text())
    assert destination.read_bytes() == data
    assert provenance["source"] == "https://example.invalid/fixture"
    assert provenance["sha256"] == digest
    assert provenance["license"] == "CC0"
    assert provenance["attribution"] == "Example Publisher"
    assert provenance["size_bytes"] == len(data)


@pytest.mark.parametrize("omitted", ["--license", "--attribution"])
def test_cli_requires_provenance_before_download(monkeypatch, tmp_path, omitted):
    def unexpected_download(url, timeout):
        pytest.fail("metadata validation must precede the download")

    monkeypatch.setattr("fieldbook.downloads.urlopen", unexpected_download)
    options = ["--sha256", "0" * 64, "--max-bytes", "100",
               "--license", "CC0", "--attribution", "Example Publisher",
               "--out", str(tmp_path / "fixture")]
    offset = options.index(omitted)
    del options[offset:offset + 2]
    with pytest.raises(SystemExit) as error:
        main(["fetch", "https://example.invalid/fixture", *options])
    assert error.value.code == 2
    assert not (tmp_path / "fixture").exists()
