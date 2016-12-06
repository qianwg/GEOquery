require(rentrez)
require(purrr)
require(lubridate)
require(xml2)

.docSumListConvert = function(x) {
    ret = data.frame(
        uid          = x %>% purrr::map('uid') %>% purrr::flatten_chr() %>% as.integer(),
        accession    = x %>% purrr::map('accession') %>% purrr::flatten_chr(),
        seriestitle  = x %>% purrr::map('seriestitle') %>% purrr::flatten_chr(),
        platformtitle= x %>% purrr::map('platformtitle') %>% purrr::flatten_chr(),
        taxon        = x %>% purrr::map('gpl') %>% purrr::map(strsplit,';') %>% purrr::flatten() %>% I(),
        entrytype    = x %>% purrr::map('entrytype') %>% purrr::flatten_chr(),
        gdstype      = x %>% purrr::map('gdstype') %>% purrr::flatten_chr(),
        ptechtype    = x %>% purrr::map('ptechtype') %>% purrr::flatten_chr(),
        valtype      = x %>% purrr::map('valtype') %>% purrr::flatten_chr(),
        ssinfo       = x %>% purrr::map('ssinfo') %>% purrr::map(strsplit,';') %>% purrr::flatten() %>% I(),
        title        = x %>% purrr::map('title') %>% purrr::flatten_chr(),
        summary      = x %>% purrr::map('summary') %>% purrr::flatten_chr(),
        gpl          = x %>% purrr::map('gpl') %>% purrr::map(strsplit,';') %>% purrr::flatten() %>% purrr::map(function(x) paste0('GPL',x)) %>% I(),
        gse          = x %>% purrr::map('gse') %>% purrr::map(strsplit,';') %>% purrr::flatten() %>% purrr::map(function(x) paste0('GSE',x)) %>% I(),
        gds          = x %>% purrr::map('gds') %>% purrr::map(strsplit,';') %>% purrr::flatten() %>% purrr::map(function(x) ifelse(nchar(x)>0,paste0('GDS',x),"")) %>% I(),
        samples      = x %>% purrr::map('samples') %>% purrr::map('accession') %>% I(),
        suppfile     = x %>% purrr::map('suppfile') %>% purrr::map(strsplit,', ') %>% purrr::flatten() %>% I(),
        ftplink      = x %>% purrr::map('ftplink') %>% purrr::flatten_chr(),
        n_samples    = x %>% purrr::map('n_samples') %>% as.integer(),
        pubmedids    = x %>% purrr::map('pubmedids') %>% I(),
        projects     = x %>% purrr::map('projects') %>% I(),
        public_date  = x %>% purrr::map('pdat') %>% purrr::flatten_chr() %>% lubridate::date(),
        samplestaxa  = x %>% purrr::map('samplestaxa') %>% I() #purrr::map(strsplit,', ') %>% purrr::flatten() %>% I()
        )
    return(ret)
    }

.docSumListFix = function(x) {
  if(!is.list(x)) x = list(x)
  ret = lapply(x,function(val) {
    val$gpl = list(paste0('GPL',stringi::stri_split(val$gpl,fixed=';')[[1]]))
    val$gse = list(paste0('GSE',stringi::stri_split(val$gse,fixed=';')[[1]]))
    val$gds = list(paste0('GDS',stringi::stri_split(val$gds,fixed=';')[[1]]))
    val$ssinfo  = unique(stringi::stri_split(val$ssinfo,fixed=';')[[1]])
    return(val)
  })
  return(ret)
}

#' get all GEO records as docsums
#'
#' Fetches all records from NCBI entrez and creates a data frame
#' from the docsum entries, loadable into a SQL-like database
#'
#'
searchGEODocSums = function(term='GPL[ETYP] OR GSE[ETYP] OR GSM[ETYP] OR GDS[ETYP]') {
    res = rentrez::entrez_search(db='gds',use_history=TRUE,term=term)
    return(res)
}
fetchDocSums= function(res,retstart=1,retmax=100) {
  z=NULL
  while(is.null(z) || inherits(z,'error')) {
    z = try(rentrez::entrez_summary(db='gds',
                                    web_history=res$web_history,
                                    retstart=retstart,
                                    retmax=retmax))
  }
    
  return(.docSumListConvert(z))
}
.bqSchemaFromDocsums = function(exampleDF) {
  mappings = list(integer=list(mode='nullable',type='integer'),
                  factor =list(mode="nullable",type="string"),
                  character = list(mode='nullable',type='string'),
                  numeric = list(mode="nullable",type="float"),
                  logical = list(mode="nullable",type="boolean"),
                  AsIs    = list(mode="repeated",type="string"),
                  Date    = list(mode="nullable",type="string"))
  coltypes = sapply(exampleDF,class)
  coldesc  = mappings[coltypes]
  names(coldesc)=NULL
  coldesc=lapply(seq_along(coldesc),function(col) {
    ret = coldesc[[col]]
    ret$name=colnames(exampleDF)[col]
    return(ret)
  })
  return(coldesc)
}


getGEOMeta = function(geo) {
    .getGSMmeta = function(dat) {
      columns = data.frame(do.call(rbind,lapply(xml_find_all(dat,'/d1:MINiML/d1:Sample/d1:Data-Table/d1:Column'),
                          function(column) {
                            return(list(
                              name = xml_text(xml_find_first(column,'d1:Name')),
                              description = xml_text(xml_find_first(column,'d1:Description'))

                            ))
                          })))
        ret = list(
          details = data.frame(
            accession     = xml_attr(xml_find_first(dat,'/d1:MINiML/d1:Sample'),'iid'),
            gpl           = xml_attr(xml_find_first(dat,'/d1:MINiML/d1:Platform'),'iid'),
            organization  = xml_text(xml_find_first(dat,'/d1:MINiML/d1:Contributor/d1:Organization')),
            title         = xml_text(xml_find_first(dat,'/d1:MINiML/d1:Sample/d1:Title')),
            description   = xml_text(xml_find_first(dat,'/d1:MINiML/d1:Sample/d1:Description'),trim = TRUE),
            type          = xml_text(xml_find_first(dat,'/d1:MINiML/d1:Sample/d1:Type')),
            submission_date= date(xml_text(xml_find_first(dat,'/d1:MINiML/d1:Sample/d1:Status/d1:Submission-Date'))),
            last_update   = date(xml_text(xml_find_first(dat,'/d1:MINiML/d1:Sample/d1:Status/d1:Last-Update-Date'))),
            release_date  = date(xml_text(xml_find_first(dat,'/d1:MINiML/d1:Sample/d1:Status/d1:Release-Date'))),
            n_channels    = as.integer(xml_text(xml_find_first(dat,'/d1:MINiML/d1:Sample/d1:Channel-Count')))
          ),
            channels      = data.frame(do.call(rbind,lapply(xml_find_all(dat,'/d1:MINiML/d1:Sample/d1:Channel'),
                                   function(channel) {
                                       return(list(
                                           tax_id = xml_attr(xml_find_first(channel,'./d1:Organism'),'taxid'),
                                           organism = xml_text(xml_find_first(channel,'./d1:Organism')),
                                           source = xml_text(xml_find_first(channel,'./d1:Source')),
                                           molecule = xml_text(xml_find_first(channel,'./d1:Molecule')),
                                           treatment_protocol = xml_text(xml_find_first(channel,'./d1:Treatment-Protocol')),
                                           extract_protocol = str_trim(xml_text(xml_find_first(channel,'./d1:Extract-Protocol'))),
                                           label_protocol = xml_text(xml_find_first(channel,'./d1:Label-Protocol')),
                                           scan_protocol = xml_text(xml_find_first(channel,'./d1:Scan-Protocol')),
                                           hybridization_protocol = xml_text(xml_find_first(channel,'./d1:Hybridization-Protocol')),
                                           growth_protocol = str_trim(xml_text(xml_find_first(channel,'./d1:Growth-Protocol'))),
                                           characteristics = sapply(xml_find_all(channel,'d1:Characteristics'),function(x) {str_trim(xml_text(x))}) %>%
                                             setNames(sapply(xml_find_all(channel,'d1:Characteristics'),function(x) {str_trim(xml_attrs(x,'tag'))})) %>% I()
                                       ))
                                   }))),
          columns = data.frame(do.call(rbind,lapply(xml_find_all(dat,'/d1:MINiML/d1:Sample/d1:Data-Table/d1:Column'),
                                                    function(column) {
                                                      return(list(
                                                        name = xml_text(xml_find_first(column,'d1:Name')),
                                                        description = xml_text(xml_find_first(column,'d1:Description'))
                                                        
                                                      ))
                                                    })))
                               
        )
        return(ret)
    }
    .getGSEmeta = function(dat) {
      pubmed = as.integer(xml_text(xml_find_all(dat,'/d1:MINiML/d1:Series/d1:Pubmed-ID')))
      pubmed = ifelse(length(pubmed)==0,NA,pubmed)
      ret = data.frame(
        accession     = xml_attr(xml_find_first(dat,'/d1:MINiML/d1:Series'),'iid'),
        samples       = list(xml_text(xml_find_all(dat,'d1:Sample'))) %>% I(),
        organization  = xml_text(xml_find_first(dat,'/d1:MINiML/d1:Contributor/d1:Organization')),
        gpl           = xml_attr(xml_find_first(dat,'/d1:MINiML/d1:Platform'),'iid'),
        title         = xml_text(xml_find_first(dat,'/d1:MINiML/d1:Series/d1:Title')),
        summary       = xml_text(xml_find_first(dat,'/d1:MINiML/d1:Series/d1:Summary')),
        overall_design= xml_text(xml_find_first(dat,'/d1:MINiML/d1:Series/d1:Overall-Design')),
        type          = xml_text(xml_find_first(dat,'/d1:MINiML/d1:Series/d1:Type')),
        pubmed        = pubmed,
        submission_date= date(xml_text(xml_find_first(dat,'/d1:MINiML/d1:Series/d1:Status/d1:Submission-Date'))),
        last_update   = date(xml_text(xml_find_first(dat,'/d1:MINiML/d1:Series/d1:Status/d1:Last-Update-Date'))),
        release_date  = date(xml_text(xml_find_first(dat,'/d1:MINiML/d1:Series/d1:Status/d1:Release-Date')))
      )
      if(!is.data.frame(ret)) return(NULL)
      return(ret)
    }
    .getGPLmeta = function(dat) {
      ret = data.frame(
        accession     = xml_attr(xml_find_first(dat,'/d1:MINiML/d1:Platform'),'iid'),
        organization  = xml_text(xml_find_first(dat,'/d1:MINiML/d1:Contributor/d1:Organization')),
        title         = xml_text(xml_find_first(dat,'/d1:MINiML/d1:Platform/d1:Title')),
        description   = xml_text(xml_find_first(dat,'/d1:MINiML/d1:Platform/d1:Description')),
        manufacturer  = xml_text(xml_find_first(dat,'/d1:MINiML/d1:Platform/d1:Manufacturer')),
        technology    = xml_text(xml_find_first(dat,'/d1:MINiML/d1:Platform/d1:Technology')),
        protocol      = xml_text(xml_find_first(dat,'/d1:MINiML/d1:Platform/d1:Manufacture-Protocol')),
        submission_date= date(xml_text(xml_find_first(dat,'/d1:MINiML/d1:Platform/d1:Status/d1:Submission-Date'))),
        last_update   = date(xml_text(xml_find_first(dat,'/d1:MINiML/d1:Platform/d1:Status/d1:Last-Update-Date'))),
        release_date  = date(xml_text(xml_find_first(dat,'/d1:MINiML/d1:Platform/d1:Status/d1:Release-Date')))
      )
      if(!is.data.frame(ret)) return(NULL)
      return(ret)
    }
    geo = toupper(geo)
    stopifnot(substr(geo,1,3) %in% c('GSE','GSM','GPL'))
              
    dat = read_xml(sprintf("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s&view=brief&targ=self&form=xml",geo))
    ret = switch(substr(geo,1,3),
                 GSM = .getGSMmeta(dat),
                 GSE = .getGSEmeta(dat),
                 GPL = .getGPLmeta(dat),
                 NA)
    return(ret)
}


