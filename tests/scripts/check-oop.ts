// Inspect the OOP parameter table: wildcard coverage per central type.
import { OOP_PARAMS } from '../../src/mmff94/parameters/out-of-plane';

const keys = Object.keys(OOP_PARAMS);
console.log('total keys:', keys.length);
// The wildcard shape: 5 parts, "cls-i-j-k-l" with i=k=l=0 and j = central
const wildcards = keys.filter(k => {
  const p = k.split('-');
  return p.length === 5 && p[1] === '0' && p[3] === '0' && p[4] === '0';
});
console.log('wildcard central types:', wildcards.map(k => k.split('-')[2]).join(','));
console.log('40-wildcard present:', wildcards.some(k => k.split('-')[2] === '40'));
