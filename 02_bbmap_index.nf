include {bbmap_index} from  "${projectDir}/subscripts/modules.nf"
import static groovy.io.FileType.FILES
import groovy.io.FileType


ref_list = params.ref.split(/,/)
c = 1
//finding if there builded indexes in working directory
def dir = new File("${workflow.launchDir}/${params.work_dir}/ref/genome")
if ( dir.exists() ) {
    def list = []
    dir.eachFileRecurse (FileType.DIRECTORIES) { file ->
        list << file.name.split(/\//)[-1]
    }
    c = list.max().toInteger() + 1
}

workflow {
    Channel
        .fromPath( ref_list.toList() )
        .map { tuple( it, c++ ) }
        .set { ref_ch }
    bbmap_index( ref_ch )
}
