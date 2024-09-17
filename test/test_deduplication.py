from pathlib import Path

import pytest

from general_utilities.job_management.command_executor import DockerMount, CommandExecutor
from makebgen.deduplication.deduplication import deduplicate_variants



@pytest.mark.parametrize(argnames=['vep_id', 'previous_vep_id', 'vcf_prefix'],
                         argvalues=[('file-GqBqggjJx5qp6yV4px0ZYgJ7', 'file-GqBqkvQJx5qj70P21QpZKf7Q', Path('ukb24310_c7_b5011_v1'))])
def test_deduplicate_variants(tmp_path, vep_id, previous_vep_id, vcf_prefix):

    test_mount = DockerMount(tmp_path, Path('/test/'))
    cmd_exec = CommandExecutor(docker_image='egardner413/mrcepid-burdentesting', docker_mounts=[test_mount])

    deduped = deduplicate_variants(vep_id, previous_vep_id, vcf_prefix, cmd_exec)
    assert deduped.exists()



