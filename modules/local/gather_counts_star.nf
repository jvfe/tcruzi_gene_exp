process GATHER_COUNTS_STAR {
    label "process_medium"

    container "biocontainers/r-tidyverse:1.2.1"

    input:
    path counts

    output:
    path "count_table_star.txt"
    path "count_table_star.rds"

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    get_counts_star.R \\
        $counts

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """

}
