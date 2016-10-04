#pragma once

#include <rdma/fabric.h>
#include <rdma/fi_domain.h>
#include <rdma/fi_endpoint.h>

extern struct fi_info* info;
extern struct fid_fabric *fabric;
extern struct fid_domain *domain;
extern struct fid_ep *ep;
extern struct fid_cntr *rxcntr;
extern struct fid_cntr *txcntr;
extern struct fid_av *av;

extern void setup_ofi_with_cntr(bool with_rma_event);
