api_key = "AIzaSyDQXJyoRCOprAYdZLK-p-tKT_8YfKVYrik"
llm_interpreted_functional_module <-
  llm_interpret_module(
    object = biotext_functional_modules,
    module_content_number_cutoff = 2,
    api_provider = "gemini",
    api_key = api_key,
    llm_model = "gemini-1.5-flash",
    embedding_model = "text-embedding-004",
    orgdb = org.Hs.eg.db,
    embedding_output_dir = "C:\\Users\\cgxjd\\Downloads"
  )
