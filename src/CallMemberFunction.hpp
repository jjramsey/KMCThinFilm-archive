#ifndef CALL_MEMBER_FUNCTION_HPP
#define CALL_MEMBER_FUNCTION_HPP

// From the C++ FAQ, http://www.parashift.com/c++-faq/macro-for-ptr-to-memfn.html
#define KMC_CALL_MEMBER_FUNCTION(object,ptrToMember)  ((object).*(ptrToMember))

#endif /* CALL_MEMBER_FUNCTION_HPP */
