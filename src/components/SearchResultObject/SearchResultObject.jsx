import React from "react";
import {Col, Container, Row, Badge, Placeholder, Popover, OverlayTrigger, Alert} from "react-bootstrap";
import './SearchResultObject.css'
import DrugDescription from "../DrugDescription/DrugDescription";
import DrugAccordion from "../DrugAccordion/DrugAccordion";
import ChartDisplayObject from "../ChartDisplayObject/ChartDisplayObject";


const SearchResultObject = (props) => {
    const { termCountDict, totalCount } = props.searchResults.result
    const eventsOverTime = props.eventsOverTime.result.results
    const [retrievedTermDataArr, setRetrievedTermDataArr] = React.useState([])

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

    return (
        <div>
            <Container className={'my-4'}>
                <Row>
                    {retrievedTermDataArr.length > 0 ? (
                            retrievedTermDataArr.map((term, index) => (
                                <Col key={index}>
                                    <h1>{term.drug_name}</h1>
                                    <h3>{processDrugGroups(term.groups)}</h3>
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
                <Row>
                    <Alert variant={"secondary"}>
                        <span>The provided information is indicative and shouldn't be used for inference.</span>
                    </Alert>
                </Row>
                <Row>
                    {props.searchTerm !== undefined && (
                        Array.isArray(props.searchTerm)
                        ? props.searchTerm.map((term, index) => (
                            <Col key={index}>
                                <DrugDescription drugName={term} onRetrieved={setRetrievedTermDataArr}
                                                 showAdditionalSearch={props.showAdditionalSearch}/>
                            </Col>
                        ))
                        : (
                            <Col>
                                <DrugDescription drugName={props.searchTerm} onRetrieved={setRetrievedTermDataArr}
                                                 showAdditionalSearch={props.showAdditionalSearch}/>
                            </Col>
                        )
                    )}
                </Row>
                <Row>
                    {Array.isArray(retrievedTermDataArr) && retrievedTermDataArr.length > 0 && (
                        retrievedTermDataArr.map((term, index) => (
                            <Col key={index}>
                                <DrugAccordion retrievedTermArr={term} totalCount={totalCount}/>
                            </Col>
                        )))
                    }
                </Row>
                <Row>
                    <Col className={'mt-3'}>
                        <ChartDisplayObject eventsOverTime={eventsOverTime} termCountDict={termCountDict}
                                            totalCount={totalCount} searchResults={termCountDict}/>
                    </Col>
                </Row>
            </Container>
        </div>
    );
}

export default SearchResultObject;