from unittest.mock import Mock

from sftool.variant_confirmation.liftover import GeneBeVariantLiftover
from sftool.variant_confirmation.models import VariantCandidate


def test_lift_grch38_candidate_to_grch37():
    session = Mock()
    response = Mock()

    response.raise_for_status.return_value = None
    response.json.return_value = {
        "dest": "HG19",
        "from": "HG38",
        "variants": [
            {
                "alt": "T",
                "chr": "chr7",
                "pos": 117199641,
                "ref": "TATC",
            }
        ],
    }

    session.get.return_value = response

    liftover = GeneBeVariantLiftover(
        session=session,
        timeout=30.0,
    )

    candidate = VariantCandidate(
        candidate_id="candidate_1",
        chromosome="7",
        position=117559587,
        reference="TATC",
        alternate="T",
        assembly="GRCh38",
        conversion_warnings=[],
    )

    lifted = liftover.lift(
        candidates=[candidate],
        target_assembly="GRCh37",
    )

    assert len(lifted) == 1

    result = lifted[0]

    assert result.candidate_id == "candidate_1"
    assert result.chromosome == "chr7"
    assert result.position == 117199641
    assert result.reference == "TATC"
    assert result.alternate == "T"
    assert result.assembly == "GRCh37"

    session.get.assert_called_once_with(
        liftover.API_URL,
        params={
            "query": "7-117559587-TATC-T",
            "from": "hg38",
            "dest": "hg19",
        },
        headers={
            "Accept": "*/*",
        },
        timeout=30.0,
    )