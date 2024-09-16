from pathlib import Path

import pytest

from makebgen.deduplication.deduplication import deduplicate_variants

@pytest.mark.parametrize(argnames=['vep_id', 'previous_vep_id', 'vcf_prefix'],
                         argvalues=[('file-GqBqggjJx5qp6yV4px0ZYgJ7', 'file-GqBqkvQJx5qj70P21QpZKf7Q', Path('ukb24310_c7_b5011_v1'))])
def test_deduplicate_variants(vep_id, previous_vep_id, vcf_prefix):

    deduped = deduplicate_variants(vep_id, previous_vep_id, vcf_prefix)
    assert deduped.exists()



