#include "ofi.hpp"

#include <rdma/fi_errno.h>
#include <string.h>

#include <cassert>
#include <cstdio>
#include <iostream>

struct fi_info* info = nullptr;
struct fid_fabric *fabric = nullptr;
struct fid_domain *domain = nullptr;
struct fid_ep *ep = nullptr;
struct fid_cntr *rxcntr = nullptr;
struct fid_cntr *txcntr = nullptr;
struct fid_av *av = nullptr;
const uint32_t num_cqe_ = 1000000;

using namespace std;

void setup_ofi_with_cntr(bool with_rma_event) {
  struct fi_info* hints = fi_allocinfo();

  if(with_rma_event)
    hints->caps = FI_RMA | FI_WRITE | FI_REMOTE_WRITE | FI_RMA_EVENT | FI_ATOMIC;
  else
    hints->caps = FI_RMA | FI_WRITE | FI_REMOTE_WRITE | FI_ATOMIC;
  hints->ep_attr->type = FI_EP_RDM;
  hints->domain_attr->control_progress = FI_PROGRESS_AUTO;
  hints->domain_attr->data_progress = FI_PROGRESS_AUTO;
  hints->domain_attr->threading = FI_THREAD_SAFE;
  hints->domain_attr->mr_mode = FI_MR_BASIC;
  hints->addr_format = FI_SOCKADDR_IN;
  hints->fabric_attr->prov_name = "gni";
  //hints->fabric_attr->prov_name = "sockets";

  // getinfo
  int res = fi_getinfo(FI_VERSION(1, 1), nullptr, nullptr, 0, hints, &info);
  assert(res == 0);

  // fabric
  res = fi_fabric(info->fabric_attr, &fabric, NULL);
  assert(res == 0);

  //domain
  res = fi_domain(fabric, info, &domain, NULL);
  assert(res == 0);

  // endpoint
  res = fi_endpoint(domain, info, &ep, NULL);
  assert(res == 0);

  // counter
  struct fi_cntr_attr tx_cntr_attr, rx_cntr_attr;
  memset(&tx_cntr_attr, 0, sizeof(tx_cntr_attr));
  tx_cntr_attr.events = FI_CNTR_EVENTS_COMP;
  tx_cntr_attr.wait_obj = FI_WAIT_UNSPEC;
  res = fi_cntr_open(domain, &tx_cntr_attr, &txcntr, NULL);
  assert(res == 0);
  res = fi_ep_bind(ep, (fid_t)txcntr, FI_READ | FI_WRITE);
  assert(res == 0);

  memset(&rx_cntr_attr, 0, sizeof(rx_cntr_attr));
  rx_cntr_attr.events = FI_CNTR_EVENTS_COMP;
  rx_cntr_attr.wait_obj = FI_WAIT_UNSPEC;
  res = fi_cntr_open(domain, &rx_cntr_attr, &rxcntr, NULL);
  assert(res == 0);
  res = fi_ep_bind(ep, (fid_t)rxcntr, FI_REMOTE_WRITE);
  assert(res == 0);

  // av
  struct fi_av_attr av_attr;
  memset(&av_attr, 0, sizeof(av_attr));
  av_attr.type = FI_AV_TABLE;
  av_attr.count = 100;
  res = fi_av_open(domain, &av_attr, &av, NULL);
  assert(res == 0);
  res = fi_ep_bind(ep, (fid *)av, 0);
  assert(res == 0);

  res = fi_enable(ep);
  if(res != 0)
    printf("fi_enable %s\n", fi_strerror(-res));
  assert(res == 0);
}
