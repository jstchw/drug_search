import React from "react";
import {Col, Container, Row, Badge, Placeholder, Popover, OverlayTrigger, Alert} from "react-bootstrap";
import './SearchResultObject.css'
import DrugDescription from "../DrugDescription/DrugDescription";
import DrugAccordion from "../DrugAccordion/DrugAccordion";
import ChartDisplayObject from "../ChartDisplayObject/ChartDisplayObject";
import PatientCard from "../PatientCard/PatientCard";
import {searchTypes} from "../OptionModal/OptionModal";
import useDrugInfo from "../../hooks/useDrugInfo";

const getUniqueObjects = (array) => {
    return array.filter((obj, index, self) =>
        index === self.findIndex((otherObj) => obj.drug_name === otherObj.drug_name)
    );
}


const SearchResultObject = (props) => {
    const { termCountDict, totalCount } = props.searchResults.result
    const eventsOverTime = props.eventsOverTime.result.results
    const searchTypeRef = React.useRef(props.searchOptions.searchBy)

    const drugInfo = useDrugInfo(props.searchTerm, searchTypeRef.current)

    const uniqueDrugInfo = getUniqueObjects(drugInfo)
    let groupDescription

    const processDrugGroups = (groups) => {
        return groups.map((group, index) => {
            let variant
            switch (group) {
                case 'approved':
                    variant = 'success'
                    groupDescription = 'Has been approved in at least one jurisdiction, at some point in time.'
                    break
                case 'investigational':
                    variant = 'warning'
                    groupDescription = 'Undergoing evaluation in the drug approval process in at least one jurisdiction.'
                    break
                case 'illicit':
                    variant = 'danger'
                    groupDescription = 'Deemed illegal in at least one jurisdiction, at some point in time.'
                    break
                case 'vet_approved':
                    variant = 'info'
                    groupDescription = 'Approved for veterinary use in at least one jurisdiction, at some point in time.'
                    break
                case 'withdrawn':
                    variant = 'secondary'
                    groupDescription = 'Has been withdrawn from the market in at least one jurisdiction, at some point in time.'
                    break
                case 'nutraceutical':
                    variant = 'light'
                    groupDescription = 'Pharmaceutical-grade nutrient with potential health benefits.'
                    break
                case 'experimental':
                    variant = 'dark'
                    groupDescription = 'Shown to bind specific proteins in experimental settings.'
                    break
                default:
                    variant = 'secondary'
            }
            group = group.substring(0, 3).toUpperCase()

            const popover = (
                <Popover id={`popover-${index}`}>
                    <Popover.Body>
                        <div className={'fs-6 text-center'}>{groupDescription}</div>
                    </Popover.Body>
                </Popover>
            )

            return (
                <OverlayTrigger key={index} trigger={['hover', 'focus']} placement={'bottom'} overlay={popover}>
                    <Badge className={'mx-1 no-pointer-badge'} key={index} bg={variant}>{group}</Badge>
                </OverlayTrigger>
            )
        })

    }

    const groupedByProduct = drugInfo.reduce((acc, term) => {
        if (!acc[term.product]) {
            acc[term.product] = [];
        }
        acc[term.product].push(term);
        return acc;
    }, {});

    return (
        <div>
            <Container className={'my-4'}>
                <Row>
                    {drugInfo.length > 0 ? (
                        Object.entries(groupedByProduct).map(([product, terms], index) => (
                            <Col key={index}>
                                {/* Searching by active substance*/}
                                {searchTypeRef.current === searchTypes[0].value &&
                                    <>
                                        <h1>{terms[0].drug_name}</h1>
                                        {terms[0].groups && <h3>{processDrugGroups(terms[0].groups)}</h3>}
                                    </>
                                }
                                {/* Searching by brand names*/}
                                {searchTypeRef.current === searchTypes[1].value &&
                                    <>
                                        <h1>{product}</h1>
                                        <Badge className={'mb-3'}>Active substance</Badge>
                                        {terms.map((term, termIndex) => (
                                            <React.Fragment key={termIndex}>
                                                <h5>{term.drug_name}{processDrugGroups(term.groups)}</h5>
                                            </React.Fragment>
                                        ))}
                                    </>
                                }
                                {/* Searching by adverse effect*/}
                                {searchTypeRef.current === searchTypes[2].value &&
                                    <>
                                        <h1>{terms.map(term => term.ADE).join(' & ')}</h1>
                                    </>
                                }
                            </Col>
                        ))
                        ) : (
                            <>
                                <Placeholder as={'h1'} animation="glow">
                                    <Placeholder xs={6} />
                                </Placeholder>
                                <Placeholder as={'h3'} animation="glow">
                                    <Placeholder xs={3} />
                                </Placeholder>
                            </>
                        )}
                </Row>
            </Container>

            <Container className={'main-drug-container'}>
                <React.Fragment>
                    <Alert variant={"secondary"}>
                        <span>The provided information is indicative and shouldn't be used for inference.</span>
                    </Alert>
                </React.Fragment>
                <Row>
                    {
                        (uniqueDrugInfo.length > 0 && uniqueDrugInfo
                            .every((term) => term.ADE === null && term.full_info === true)) &&
                        (uniqueDrugInfo.map((term, index) => (
                            uniqueDrugInfo.length > 1 ? (
                                <Col key={index}>
                                    <DrugDescription drugInfo={term} showAdditionalSearch={props.showAdditionalSearch}/>
                                    <DrugAccordion drugInfo={term} totalCount={totalCount}/>
                                </Col>
                            ) : (
                                <React.Fragment key={index}>
                                    <Col>
                                        <DrugDescription drugInfo={term} showAdditionalSearch={props.showAdditionalSearch}/>
                                    </Col>
                                    <Col>
                                        {/* FIX THIS ACCORDION WHEN THERE IS NO INFORMATION - MAKE IT NOT RENDER IF THERE IS NO INFO */}
                                        {/* RENDER SOME PARTS OF THE ACCORDION IF THERE IS LIMITED INFORMATION */}
                                        <DrugAccordion drugInfo={term} totalCount={totalCount}/>
                                    </Col>
                                </React.Fragment>
                            )
                        ))
                    )}
                </Row>

                <PatientCard searchOptions={props.searchOptions} totalADE={totalCount} />

            </Container>

            <Col className={'mt-3'}>
                <ChartDisplayObject eventsOverTime={eventsOverTime} termCountDict={termCountDict}
                                    totalCount={totalCount} searchResults={termCountDict}
                                    searchOptions={props.searchOptions}
                />
            </Col>
        </div>
    );
}
export default SearchResultObject;