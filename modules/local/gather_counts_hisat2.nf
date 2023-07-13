process GATHER_COUNTS_HISAT2 {
    label "process_medium"

    container "biocontainers/r-tidyverse:1.2.1"

    input:
    path feature_counts

    output:
    path "count_table_hisat2.txt"
    path "count_table_hisat2.rds"

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    get_counts_hisat2.R \\
        $feature_counts

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """

}
