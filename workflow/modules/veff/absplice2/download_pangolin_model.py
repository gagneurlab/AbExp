# Downloads one published Pangolin weight file with abexp.pangolin.download_model. The package pins
# the commit of the Pangolin repository and checks the SHA-256 sum of the file.
import logging

from abexp.pangolin import download_model

# force: an import of abexp.pangolin has already given the root logger a handler for stderr.
# The root logger stays at WARNING, so the INFO lines of polars-bio stay out of the log.
logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.WARNING,
    format="%(asctime)s %(name)s %(levelname)s: %(message)s",
    force=True,
)
logging.getLogger("abexp").setLevel(logging.INFO)
download_model(snakemake.wildcards["pangolin_model"], snakemake.params["models_dir"])
