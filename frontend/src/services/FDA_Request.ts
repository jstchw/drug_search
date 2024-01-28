// const demographicDataCache: {[key: string]: any} = {}

// export const getSideEffectsForDemographics = async (drug: string) => {
//     const demographicGroups: DemographicGroups[] = [
//         { name: 'Young Boys', sex: '1', age: [0, 18], def: 'Males aged 0-18 years'},
//         { name: 'Young Girls', sex: '2', age: [0, 18], def: 'Females aged 0-18 years'},
//         { name: 'Men', sex: '1', age: [19, 59], def: 'Males aged 19-59 years'},
//         { name: 'Women', sex: '2', age: [19, 59], def: 'Females aged 19-59 years'},
//         { name: 'Elderly Men', sex: '1', age: [60, 120], def: 'Males aged 60-120 years'},
//         { name: 'Elderly Women', sex: '2', age: [60, 120], def: 'Females aged 60-120 years'},
//     ];
//
//     // If the data for this drug is already in the cache, return it
//     if (demographicDataCache[drug]) {
//         return demographicDataCache[drug];
//     }
//
//     const results: {[key: string]: any} = {};
//
//     for (const group of demographicGroups) {
//         const searchOptions = {
//             searchBy: {
//                 value: 'patient.drug.medicinalproduct',
//                 index: 9999,
//                 label: '',
//                 type: 'searchBy',
//                 enabled: true,
//                 param: 'medicinal_product'
//             },
//             sex: {
//                 value: `patient.patientsex:${group.sex}`,
//                 index: 0,
//                 label: '',
//                 type: 'sex',
//                 enabled: true
//             },
//             age: {
//                 min: {
//                     value: group.age[0].toString(),
//                     index: 0,
//                     label: '',
//                     type: 'age'
//                 },
//                 max: {
//                     value: group.age[1].toString(),
//                     index: 1,
//                     label: '',
//                     type: 'age'
//                 },
//                 enabled: true
//             },
//             country: {
//                 value: 'US',
//                 index: 0,
//                 label: 'United States',
//                 type: 'country',
//                 enabled: false
//             }
//         };
//
//         const url = generatePath([drug], searchOptions, 'reaction')
//         const response = await fetch(url)
//         const data = await response.json()
//         if(!response.ok) {
//             break
//         }
//         results[group.name] = {
//             ...processDrugEvents(data),
//             def: group.def,
//             age: group.age
//         }
//     }
//
//     // Save the results to the cache
//     demographicDataCache[drug] = results
//
//     return results;
// };
